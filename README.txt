========================================================================
    Author: Danyang Su 
    Email: danyang.su@gmail.com
========================================================================

The code is written in Fortran and can be compiled on any platform that support Fortran. A brief introduction is shown below. For detail, please go to doc/user_guide.docx for detail. 

1. Introduction
Database conversion program for MIN3P is an efficient tool to convert geochemistry databases. Current version (V1.0.26) can convert TOUGHREACT database, CRUNCHFLOW database and PHREEQC database to MIN3P database. It can also be used to change the master variable for reactions, switch specified components and output specified data with pre-defined order. Please NOTE that NOT all the database entries are ported. Some of the parameters required by MIN3P database may be not available in other databases. In this case, these parameters are left as zero or null. Please check the converted database before using it.

2. Input files
Input files include input parameter file (*.inp), database files and alias file. The input parameter file and database files are indispensable while the alias file is optional.


3. Output files
Output files include the log file, the converted database file and name truncation data file.

4. Installation
The program is a single executable file (e.g., dbs_conversion.exe on Windows box). Double click the executable file or run through command line to start the program. Type in the input parameter file (*.inp) to start database conversion.

4.1. Run the main program dbs_conversion.exe, type in the input parameter file, e.g., “phreeqc2min3p.inp”.
4.2. The program will read in the input parameters and alias database (if exists) and then convert the database, write the output files. You may find warnings /errors during the first run.
4.3. Check the name truncation data. When the program find the name needs truncation, it will use the first 12 (can be specified in source codes) characters as the truncate name and output the original name and the truncate name into name_truncation.txt file.  You can see from the log files if there is duplicate truncate name or name conflict with the original name.  Modify the truncate name to make the name readable and then copy all the name pairs into alias.dbs.
4.4. Rerun the main program, repeat the steps 1 to 3 until all the errors are resolved. When no error occurs, there will be no name–truncation pair in the file name_truncation.txt. 
IMPORTANT: Be sure there is no error in the log file unless the error can be ignored. If there is error, try to fix it. The error information in the log file will help you.



