
**MOCCA-BH-Forge**

**Introduction:**

MOCCA-Forge is a simulation data pipeline designed to work with MOnte Carlo Cluster simulAtor (MOCCA) outputs. Currently, the script can analyze, process and dump desired data outputs using snapshot and history file outputs from MOCCA. MOCCA-Forge takes 
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
You can see the packages required to run the pipeline in simulation-pipeline/requirements.txt and install them through pip , if you are not in a conda environment.

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

   ![image](https://github.com/user-attachments/assets/515b6d03-4358-4794-a7e8-d4fa15ff99e5)

   ![image](https://github.com/user-attachments/assets/ba1a920e-f25a-4a36-915b-f211f4fdbcff)


   ![image](https://github.com/user-attachments/assets/33e8e6dd-a01f-4361-8b2d-8c9b58d434e5)


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

   **2. Cumulative counts dataframe:** After a snapshot is imported, a cumulative counts dataframe is displayed which summarizes the count of each stellar type present in the simulation at that time of evolution. For example, a snapshot of 12 Gyrs for a globular cluster is shown below:

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

   **3. Plotting the data:**A list of columns is displayed that can be used for plotting. Additionally, x, y and z projections are calculated for plotting stellar distrbution in cartesian 2D coordinates. These columns are saved as posx,posy and posz in the
   the list of available columns. Following is the brief introduction and output of some important plotting module.

   **Scatterplot Module**

   You can choose to create scatterplot of your distribution with additional parameters such as color columns and stellar type constraints.

   
   
