#! /bin/bash

echo 'Available choices:'
cd input/
ls -1
read -p 'What dir? ' -e FOLDER
cd ..

cp input/$FOLDER/*.ini input/$FOLDER/*.L2D ./ &>/dev/null
for FILE in $(find *.ini)
do
    ./engine $FILE
done

eog *.bmp

rm *.ini *.bmp *.L2D &>/dev/null
