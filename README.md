
**MOCCA-BH-Forge**

**Introduction:**

MOCCA-BH-Forge is a simulation data pipeline designed to work with MOnte Carlo Cluster simulAtor (MOCCA) outputs. Currently, the pipeline can analyze, process and dump desired data outputs using snapshot and history file outputs from MOCCA. MOCCA-BH-Forge takes 
history or snapshot files and enables the user to plot and save data subsets of their choice.

**Installation**

1. The pipeline can be cloned from github repository using:

```bash
git clone https://github.com/nebula-navigator/simulation-pipeline.git
```
2. Make sure you are in the directory:

```bash

cd simulation-pipeline
```

3. If you are using a conda environment, the required packages can be installed with:

```bash

conda install --file requirements.txt

```
If you are not in a conda environment, use:

``bash

pip install -r requirements.txt

```

**Quickstart Guide with Examples**

1. Once the requirements are installed, the pipeline can be run with:

```bash

python -m main.main

```
2. The pipeline is designed to work with either history files or snapshot files from any MOCCA simulation. The data files are processed with strict parsing methods so if the data file format is modified in any way, the code is bound to crash. The pipeline is divided into
   a main function and modules. Below is the detailed description of outputs produced by the pipeline while

**Processing History Files**

   Upon starting the pipeline the user is presented with self explanatory choices, if they want to work with snapshot file or a history file from MOCCA. Loading a history file generates a set of figures concerning the different merger stats of the object
   for which the history file has been generated.
 ![Figure1](https://github.com/user-attachments/assets/b93bf03e-c322-4a1c-bbe7-829647631e56)

 ![Figure2](https://github.com/user-attachments/assets/eac9d2c6-bf57-4595-aad0-d46b58590d0f)

![Figure3](https://github.com/user-attachments/assets/99b327a9-f8a9-4828-8f1f-6fe3566f1467)



   The plots with dark grid are made interactive so whenever a data point is clicked on one of them, all the necessary information about the event is displayed on the terminal. An example of such an interactive point showing
   a blackhole-blackhole collision is shown below:

   ```bash

   Event Information:
            --------------------
            EventType: COLLISION
            EncId: 17420626
            EncCase: 0
            Time (Myr): 410.731037
            Summary: *1BH,2BH=>1BH|O1,B0,S1,EX0,D0,MT1
            Id6: 1893666
            IdComp: 2802
            IdCompNew: 0
            IM: 2099481
            MassOld (Msun): 14943.0
            MassNew (Msun): 14965.0
            MassChange: 22.0
            aOld (Rsun): 0.0
            aNew (Rsun): 0.0
            eOld: 0.0
            eNew: 0.0
            starType: 14
            compType: 14
            binType: 0
   ```
   
   Each run of a history file will automatically generate gravitational wave candidates data involving IMBH, Neutron Stars and White Dwarfs. The data contains all the regular columns of history file for the selected compact binary starting from the 
   time it was formed until the component inspirals into IMBH. Each compact pair data is a saved in a separate dat file for further processing. The module also outputs mergersum.dat which is simply the counts table seen in figure 2.

   *Note: The current version of the module is only able to parse history outputs from latest version of MOCCA. To generate gw data from the old archives, please use histarchive.py script in the parent directory*

   **Processing Snapshot Files**

   **1. Working with Multiple Snapshots:** The pipeline is capable of analyzing snapshot files from MOCCA in detail. The analysis includes various plot functions with desired set of parameters that can be given for each visualization type selected. After choosing to work with snapshot file,
   you can load multiple snapshots at once and then choose to work with the one you are interested in. Due to the large size of MOCCA snapshot files > 700MB, it is not feasible to plot and visualize the data from multiple snapshot files all at once. But
   after the pipeline has loaded the files it is fairly easy to switch between them while staying in the program. You can load multiple snapshots using a comma.

   **2. Cumulative counts dataframe:** After a snapshot is imported, a cumulative counts dataframe is displayed which summarizes the count of each stellar type present in the simulation at that time of evolution. For example, cumulative counts dataframe of a snapshot of 12 Gyrs for a globular cluster is shown below:

```bash
    Cumulative Counts DataFrame:
                                          Stellar Type  Cumulative Count
    0                      Low-mass main-sequence star           1489847
    1                               Main-sequence star             59112
    2                             Hertzsprung-gap star              1011
    3                          First giant branch star              2287
    4                         Core helium burning star               696
    5               Early asymptotic giant branch star                27
    6   Thermally pulsing asymptotic giant branch star                 2
    7                             Naked helium star MS                 1
    8                Naked helium star Hertzsprung gap                16
    9                   Naked helium star giant branch                 0
    10                              Helium white dwarf               945
    11                       Carbon-oxygen white dwarf            185267
    12                         Oxygen-neon white dwarf              3470
    13                                    Neutron star              1781
    14                                      Black hole                 2
    15                              Massless supernova                 0

```

 **3. Plotting the data:** A list of columns is displayed that can be used for plotting. Additionally, x, y and z projections are calculated for plotting stellar distrbution in cartesian 2D coordinates. These columns are saved as posx,posy and 
      posz in the the list of available columns. 

   **Current Visualization Options**

   ```bash
   Suggested visualizations:
   1. Histogram: For visualizing the distribution of a single continuous variable.
   2. Scatter Plot: For visualizing the relationship between two continuous variables, grouped by 2 color axes
   3. 3D Scatter Plot: For visualizing the relationship between three continuous variables.
   4. Box Plot: For visualizing the distribution of a continuous variable, possibly grouped by a categorical variable.
   5. Pair Plot: For visualizing the relationships between multiple continuous variables.
   6. Violin Plot: For visualizing the distribution of a continuous variable, possibly grouped by a categorical variable.
   7. Heatmap: For visualizing the correlation matrix between multiple continuous variables.
   8. 3D Visualization of Positions with Interactive Plot: For visualizing the positions of objects in the cluster in 3D space with interactive features.
   9. Bar Plot of Event Frequencies: For visualizing the frequency of different events (takes single column)
   10. Cumulative Events vs Stellar Types ( bar plot of cumulative sum of each stellar type)
   11. Position vs Stellar Type (radial distribution of each stellar type
   12. Stellar type distribution by mass and radial position. A powerful plot for multiple distribution analysis
   13. Cumulative sum of each stellar type vs radial position line plot (suggested plot from reference thesis)
   14. Plot HR diagram
   15. Plot Blue Straggler cumulative counts vs r
   16. Return to Main Menu

```
   **Example plots from snapshot files**

   You can create any customized plot using visualization options. Here are a few examples of plots that can be generated with MOCCA-BH-Forge.

   ![3d-positions](https://github.com/user-attachments/assets/a698ef58-2e82-4340-89ad-0d4ec4fd2183)

   ![Figure_4](https://github.com/user-attachments/assets/1cc27a05-b5a6-4078-b5ea-7c2c379a9220)

   ![Figure5](https://github.com/user-attachments/assets/466c99a8-c3ba-4ba1-a519-ea4c0869d815)

   ![Figure6,png](https://github.com/user-attachments/assets/d8311f6c-724f-4de1-a3bc-10229740de00)




   

   
   
