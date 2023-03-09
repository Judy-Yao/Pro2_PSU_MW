#!/bin/bash
#With the input order directory (such as 8290131495), this script first download the website with file names to be downloaded, then grep file address from the website, last serially download each file.

# name of the order
order=$1

# download parent address
parent_address=https://download.avl.class.noaa.gov/download/${order}/001

# download the website 
wget $parent_address

# grep OR-related records
grep "OR" 001 > file_address

# read the file line by line
while IFS= read -r line
do
 # get the file address on that record
 file_name=$(echo "$line" |cut -c54-129)
 echo 'Downloading' $file_name
 wget ${parent_address}/${file_name}

done < file_address

