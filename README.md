
# Order of Execution (Files Only)

1. `ProteoCorrAnalyser.py`
2. `BayesianWeightOptimization.py`
3. `ProximityBasedCysFeSelector.py`
4. `ANMModeAveragesCompiler.py`
5. `island_dynamic_programming_V9.py`
6. `TileBasedSimilarityVisualization.R`
7. `ProteinTileSimViz.R`
8. `ProteinChainOptimalPathAnalyzer.R`

## Example execution:
1\. **`ProteoCorrAnalyser.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ProteoCorrAnalyser.py 5t5i 2fdn All G A;
   ```

   #### Outputs:
   - **Cross-Correlation Data CSVs**: Files named `"{pdb_id}_CrossCorr_Mode{anm_mode_number}.csv"`. These CSV files contain cross-correlation data between protein chains, tagged with the PDB ID and the ANM mode number. Useful for analyzing the correlation in movements between different parts of a protein or between different proteins.

   - **Filtered Average Cosine Similarities Between Chains CSVs**: Files named `"filtered_average_cosine_similarities_Mode{anm_mode_number}_between_chain_{chain_id_one}_and_{chain_id_two}.csv"`. These files provide the filtered average cosine similarities between specified chains, offering insights into the similarity of movements or orientations between different chains within the same or different proteins.

   - **Cosine Similarities Within Protein CSVs**: Files named `"cosine_similarities_Mode{modeNum}_within_protein.csv"` and `"cosine_similarities_Mode{modeNum}.csv"`. These CSV files detail the cosine similarities within a protein or more broadly, helping in the assessment of how different parts of a protein or different proteins move in relation to each other.



2\. **`BayesianWeightOptimization.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/BayesianWeightOptimization.py All 5t5i.chainG.nowa
   ```

   #### `BayesianWeightOptimization.py` Outputs:

   - **Cosine Similarities Updated Ranked Local Alignment CSVs**:
      - `"cosine_similarities_Mode{mode_num}_updated_ranked_LocalAln.csv"`: Contains ranked cosine similarities after local alignment updates, providing insights into the alignment optimization process.

   - **Cosine Similarities Updated Local Alignment CSVs**:
      - `"cosine_similarities_Mode{mode_num}_updated_LocalAln.csv"`: Details the cosine similarities after updates from local alignment, highlighting adjustments made to improve alignment accuracy.

   - **Cosine Similarities Updated Ranked Local Alignment Weighted Sum CSVs**:
      - `"cosine_similarities_Mode{mode_num}_updated_ranked_LocalAln_weightedSum_whole.csv"`: Features the weighted sum of ranked cosine similarities post-local alignment update, offering a comprehensive view of alignment optimization across the entire dataset.

   - **Generic Cosine Similarities CSVs** (possibly input or intermediary files):
      - `"cosine_similarities_Mode{mode_num}_within_protein.csv"` and `"cosine_similarities_Mode{mode_num}.csv"`: These names suggest the files are related to cosine similarities, potentially serving as input or intermediary files in the optimization process.


3\. **`ProximityBasedCysFeSelector.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ProximityBasedCysFeSelector.py All 5t5i G;
   ```

   #### `ProximityBasedCysFeSelector.py` Outputs:

   - **Final Scan Mode CSVs**:
      - `"final_scan_Mode{mode_num}.csv"`: This CSV file contains the results of a final scan, likely involving the selection of cysteine and iron (Cys-Fe) pairs based on proximity or other criteria, ready for further analysis.

   - **Final Scan Mode with False CSVs**:
      - `"final_scan_Mode{mode_num}_withFalse.csv"`: Similar to the above, this file includes the results of the final scan but potentially includes additional data or flags (e.g., false positives) for detailed examination.




4\. **`ANMModeAveragesCompiler.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ANMModeAveragesCompiler.py All
   ```

   #### `ANMModeAveragesCompiler.py` Outputs:

   - **Averages and Counts Mode CSV**:
      - `"averages_and_counts_Mode{anm_mode_number}.csv"`: This CSV file contains the averages and counts data, likely related to Anisotropic Network Model (ANM) mode analysis. This data could include statistical averages and occurrence counts of certain features or metrics across the analyzed modes.

      - **Note**: The script also references reading from a CSV file named `"final_scan_Mode{anm_mode_number}_withFalse.csv"`, which might be an input file from a previous analysis step, rather than an output of this particular script.





5\. **`island_dynamic_programming_V9.py`**

   ```bash
   python3.11 /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/island_dynamic_programming_V9.py
   ```

   #### `island_dynamic_programming_new.py` Outputs:

   - **Path Sequences CSV**:
      - `"path_sequences.csv"`: This CSV file contains the sequences of paths determined by the dynamic programming algorithm, likely used for identifying optimal paths or sequences within a dataset or analysis context.

   - **All Labeled Points CSV**:
      - `"all_labeled_points.csv"`: This file includes all points (or data entries) that have been labeled through the analysis process, potentially indicating significance or categorization based on the algorithm's criteria.

   - **Note**: The script also references reading from CSV files, such as `"final_scan_ModeAll_withFalse.csv"` and `"averages_and_counts_ModeAll.csv"`, which might be input files from previous analysis steps, rather than outputs of this particular script.





6\. **`TileBasedSimilarityVisualization.R`**

   ```bash
   Rscript --vanilla /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/TileBasedSimilarityVisualization.R 5t5i "Chain G" All
   ```

7\. **`ProteinTileSimViz.R`**

   ```bash
   Rscript --vanilla /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ProteinTileSimViz.R 5t5i "Chain G" All
   ```

8\. **`ProteinChainOptimalPathAnalyzer.R`**

   ```bash
   Rscript --vanilla /Users/jansiess/Desktop/main_repository/v/vikLab/projects/thesis/aim1/paper_1_code/current_version/ProteinChainOptimalPathAnalyzer.R --pdb-id=5t5i --chain-id=G
   ```

