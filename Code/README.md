- **racipe_run.py** - Simulates each topology file in triplicates using multithreaded RACIPE and stores the output in separate directories created for each RACIPE run.
- **gknorm_multi.py** - Removes extra solution files generated by RACIPE. Names and g/k normalizes the solution, stores the output in the directories already created for each run. 
- **Correlation_multi.py** - Generates correlation matrices from the g/k normalized solution and stores them in the directory for each RACIPE run. Analyses the matrices and creates a csv file of correlation team strengths (CTS) of all networks, along with a text file of team compositions of all networks. 
- **Influence_multi.py** - Generates influence matrices from the topology file and stores them in the directory for each RACIPE run. Analyses the matrices and creates a csv file of influence team strength (ITS) of all networks, along with a text file of team compositions of all networks.
- **PCA.py** - Analyses g/k normalized solution files to create csvs of PC1 variance, PC2 variance, Number of PCs to cover 90% variance, cumulative PCA coefficients, angles of cumulative PCA coefficients, and angles of PC1, PC2 coefficients, for all networks. Also generates a text file of team compositions for all networks based on PC1,PC2 coefficient angles.  
- **Teams_HD.py** - Creates a csv file of team composition hamming distance between PCA and influence teams for all networks. 