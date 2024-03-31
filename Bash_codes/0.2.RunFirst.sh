#!/bin/bash

# The command sed -i -- 's/A/B/g' * changes all the A strings in fils in the directory to B
# This enables easy change of all the dirs in one code run

# Here, change only the directory

cd Folder_name
sed -i -- 's/Previous_Folder_name/Folder_name/g' *.txt
sed -i -- 's/Previous_Folder_name/Folder_name/g' *.txt
sed -i -- 's/Previous_Folder_name/Folder_name/g' *.sh
