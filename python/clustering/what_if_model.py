import numpy as np
import pandas as pd

from .disjoint_clusters_model import DisjointClustersModel
import util


class WhatIfModel(DisjointClustersModel):
    
    def __init__(self, max_dimensions, query_frequencies, table_name, table_scans, table_sizes, distinct_values, target_chunksize, correlations, joins, sorted_columns_during_creation, model_dir='cost_model_output/models/', model_type='boost'):
        super().__init__(max_dimensions, query_frequencies, table_name, table_scans, table_sizes, distinct_values, target_chunksize, correlations, joins, sorted_columns_during_creation)
        self.models = util.load_models(model_dir)
        self.model_formats = util.load_model_input_formats(model_dir)
        self.model_type = model_type
        self.table_scans = self.table_scans.copy()
        self.joins = self.joins.copy()
    
    
    def is_range_predicate(self, predicate_string):
        return predicate_string == "LessThan" or predicate_string == "LessThanEquals" or predicate_string == "GreaterThan" or predicate_string == "GreaterThanEquals" or "Between" in predicate_string

    def adapt_scans_to_clustering(self, scans, clustering_columns, sorting_column, dimension_cardinalities):
        # Set sortedness information
        scans['INPUT_COLUMN_SORTED'] = scans.apply(lambda x: "Ascending" if x['COLUMN_NAME'] == sorting_column else "No", axis=1)


        # Pruning - Future work: consider the following aspects
        # - first scan's input size is further reduced by the other scans
        # - if >= 2 scans: first scan will be on reference segments henceforth, second scan will be the type of the first scan


        # Guess the expected ratio of chunks with early out or all match-shortcut
        # Assumption: early out does not happen because in most cases, chunks with no matches get pruned before the shortcut could apply
        scans['NONE_MATCH_RATIO'] = 0

        # Assumption: for range queries on clustered columns, pruning yields a range where all chunks (except the two bounding chunks) match completely
        def compute_all_match_ratio(row):
            def clustering_columns_correlated_to(column):
                return [clustering_column for clustering_column in clustering_columns if column in self.correlations.get(clustering_column, {})]

            def correlates_to_clustering_column(column):
                return len(clustering_columns_correlated_to(column)) > 0

            if row['COLUMN_NAME'] in clustering_columns:
                return 1
            elif correlates_to_clustering_column(row['COLUMN_NAME']):
                return 1
            else:
                return 0

        scans['ALL_MATCH_RATIO'] = scans.apply(compute_all_match_ratio, axis=1)

        # Pruning
        scans_per_query = scans.sort_values(['INPUT_ROWS'], ascending=False).groupby(['QUERY_HASH', 'GET_TABLE_HASH'])
        for _, query_scans in scans_per_query:

            unprunable_parts = query_scans.apply(self.compute_unprunable_parts, axis=1, args=(clustering_columns, dimension_cardinalities,))
            unprunable_part = unprunable_parts.product()
            assert unprunable_part > 0, "no unprunable part"

            # Pruning: Set input size of first scan to reflect pruned chunks
            estimated_pruned_table_size = min(self.table_size, self.round_up_to_next_multiple(unprunable_part * self.table_size, self.target_chunksize))
            scans.loc[query_scans.iloc[0].name, 'INPUT_ROWS'] = np.int64(estimated_pruned_table_size)
            scans.loc[query_scans.iloc[0].name, 'INPUT_CHUNKS'] = np.int64(estimated_pruned_table_size / self.target_chunksize)

            # Pruning: if there are multiple scans and the first occurs on a clustering column, the scan is likely moved up in the PQP (i.e., executed later)
            # Consequence: Currently first scan will be performed on a reference segment, currently second scan on currently first scan's segment type
            first_scan_column_name = scans.loc[query_scans.iloc[0].name, 'COLUMN_NAME']
            if len(query_scans) > 1 and first_scan_column_name in clustering_columns:
                first_scan_segment_type = scans.loc[query_scans.iloc[0].name, 'COLUMN_TYPE']
                scans.loc[query_scans.iloc[1].name, 'COLUMN_TYPE'] = first_scan_segment_type
                scans.loc[query_scans.iloc[0].name, 'COLUMN_TYPE'] = "REFERENCE"

                # The originally first scan will be executed last, so its input size will be reduced by the other scans
                for i in range(1, len(query_scans)):
                    scans.loc[query_scans.iloc[0].name, 'INPUT_ROWS'] *= scans.loc[query_scans.iloc[i].name, 'SELECTIVITY_LEFT']

        return scans


    def estimate_table_scan_runtime(self, clustering_columns, sorting_column, dimension_cardinalities):        
        runtime = 0

        scans = self.table_scans.copy()
        scans = self.adapt_scans_to_clustering(scans, clustering_columns, sorting_column, dimension_cardinalities)
        scans = scans.drop(columns=['COLUMN_NAME', 'DESCRIPTION', 'GET_TABLE_HASH', 'LEFT_INPUT_OPERATOR_HASH', 'OPERATOR_HASH', 'OPERATOR_TYPE', 'OUTPUT_CHUNK_COUNT', 'QUERY_HASH', 'PREDICATE', 'RIGHT_INPUT_OPERATOR_HASH', 'RUNTIME_NS', 'SCANS_SKIPPED', 'SCANS_SORTED', 'TABLE_NAME', 'benefits_from_sorting', 'part_of_or_chain', 'time_per_input_row', 'time_per_output_row', 'time_per_row', 'useful_for_pruning'])
        
        scans_by_implementation = scans.groupby(['OPERATOR_IMPLEMENTATION'])
        for operator_implementation, df in scans_by_implementation:
            print(f"## Estimating {operator_implementation} scans")
            model_name = util.get_table_scan_model_name(self.model_type, operator_implementation)
            model = self.models[model_name]
            
            df = df.copy()
            df = df.drop(columns=['OPERATOR_IMPLEMENTATION'])
            df = util.preprocess_data(df)
            df = util.append_to_input_format(df, self.model_formats[model_name])
            
            predictions = model.predict(df)
            self.scan_estimates.loc[df.index, 'RUNTIME_ESTIMATE'] = np.array(predictions, dtype=np.int64)
            
            runtime += predictions.sum()
        print()
            
        return runtime

    def adapt_joins_to_clustering(self, joins, clustering_columns, sorting_column, dimension_cardinalities):
        # Set sortedness information
        def set_sortedness(row, side):
            if self.table_name != row[f"{side}_TABLE"]:
                return str(row[f"{side}_SORTED"])
            else:
                return "1" if int(row[f"{side}_SORTED"]) == 1 and row[f"{side}_COLUMN"] == sorting_column else "0"

        joins['BUILD_SORTED'] = joins.apply(set_sortedness, args=("BUILD", ), axis=1)
        joins['PROBE_SORTED'] = joins.apply(set_sortedness, args=("PROBE", ), axis=1)

        # Estimate ratio of chunks pruned by the respective GetTable operators
        joins['BUILD_TABLE_PRUNED_CHUNK_RATIO'] = 0
        joins['PROBE_TABLE_PRUNED_CHUNK_RATIO'] = 0
        scans_per_query = self.table_scans.groupby(['QUERY_HASH', 'GET_TABLE_HASH'])
        for (query_hash, _), query_scans in scans_per_query:
            unprunable_parts = query_scans.apply(self.compute_unprunable_parts, axis=1, args=(clustering_columns, dimension_cardinalities,))
            unprunable_part = unprunable_parts.product()
            assert unprunable_part >= 0 and unprunable_part <= 1, "unprunable part: " + str(unprunable_part)

            build_side_joins = joins['BUILD_TABLE'] == self.table_name
            probe_side_joins = joins['PROBE_TABLE'] == self.table_name
            query_joins = joins['QUERY_HASH'] == query_hash

            joins.loc[query_joins & build_side_joins, 'BUILD_TABLE_PRUNED_CHUNK_RATIO'] = 1 - unprunable_part
            joins.loc[query_joins & probe_side_joins, 'PROBE_TABLE_PRUNED_CHUNK_RATIO'] = 1 - unprunable_part

            # Apply pruning to joins (in some cases, scans are executed above semi joins
            # TODO: correlations?
            if unprunable_part < 1:
                clustered_scans = query_scans.apply(lambda x: x['COLUMN_NAME'] in clustering_columns, axis=1)
                clustered_scan_operator_ids = query_scans[clustered_scans]['OPERATOR_ID'].unique()

                min_clustered_scan_id = min(clustered_scan_operator_ids) if clustered_scans.any() else 0
                joins_before_clustered_scans = joins['OPERATOR_ID'] < min_clustered_scan_id

                joins.loc[query_joins & build_side_joins & joins_before_clustered_scans, 'BUILD_TABLE_ROW_COUNT'] *= unprunable_part
                joins.loc[query_joins & probe_side_joins & joins_before_clustered_scans, 'PROBE_TABLE_ROW_COUNT'] *= unprunable_part

        return joins
        
    def estimate_join_runtime(self, clustering_columns, sorting_column, dimension_cardinalities):
        runtime = 0
        
        joins = self.joins.copy()
        joins = self.adapt_joins_to_clustering(joins, clustering_columns, sorting_column, dimension_cardinalities)
        
        joins_by_implementation = joins.groupby(['OPERATOR_IMPLEMENTATION', 'BUILD_COLUMN_TYPE', 'PROBE_COLUMN_TYPE'])
        for (implementation, build_type, probe_type), df in joins_by_implementation:
            print(f"## Estimating {build_type} {probe_type} joins")
            model_name = util.get_join_model_name(self.model_type, implementation, build_type, probe_type)
            model = self.models[model_name]
            
            df = df.copy()
            df = df.drop(columns=['OPERATOR_IMPLEMENTATION', 'PREDICATE_COUNT', 'PRIMARY_PREDICATE'])            
            df = util.preprocess_data(df)
            df = df.rename(columns={
                'PROBE_TABLE_ROW_COUNT': 'PROBE_INPUT_ROWS',
                'BUILD_TABLE_ROW_COUNT': 'BUILD_INPUT_ROWS',
                'OUTPUT_ROW_COUNT': 'OUTPUT_ROWS',
                'BUILD_SORTED_0': 'BUILD_INPUT_COLUMN_SORTED_No',
                'BUILD_SORTED_1': 'BUILD_INPUT_COLUMN_SORTED_Ascending',
                'PROBE_SORTED_0': 'PROBE_INPUT_COLUMN_SORTED_No',
                'PROBE_SORTED_1': 'PROBE_INPUT_COLUMN_SORTED_Ascending'
            })

            df = df.drop(columns=['BUILDING_NS', 'BUILD_COLUMN', 'BUILD_SIDE', 'BUILD_SIDE_MATERIALIZING_NS', 'BUILD_TABLE', 'CLUSTERING_NS', 'DESCRIPTION', 'IS_FLIPPED', 'LEFT_COLUMN_NAME', 'LEFT_COLUMN_TYPE_DATA', 'LEFT_COLUMN_TYPE_REFERENCE', 'LEFT_INPUT_OPERATOR_HASH', 'LEFT_TABLE_CHUNK_COUNT', 'LEFT_TABLE_NAME', 'LEFT_TABLE_ROW_COUNT', 'OPERATOR_HASH', 'OPERATOR_TYPE', 'OUTPUT_CHUNK_COUNT', 'OUTPUT_WRITING_NS', 'PROBE_COLUMN', 'PROBE_SIDE', 'PROBE_SIDE_MATERIALIZING_NS', 'PROBE_TABLE', 'PROBING_NS', 'QUERY_HASH', 'RADIX_BITS', 'RIGHT_COLUMN_NAME', 'RIGHT_COLUMN_TYPE_REFERENCE', 'RIGHT_INPUT_OPERATOR_HASH', 'RIGHT_TABLE_CHUNK_COUNT', 'RIGHT_TABLE_NAME', 'RIGHT_TABLE_ROW_COUNT', 'RUNTIME_NS'], errors="ignore")
            df = util.append_to_input_format(df, self.model_formats[model_name])
            
            predictions = model.predict(df)
            self.join_estimates.loc[df.index, 'RUNTIME_ESTIMATE'] = np.array(predictions, dtype=np.int64)

            runtime += predictions.sum()
        print()
                        
        return runtime
    