#! /bin/bash

echo 'Available choices:'
cd input/
ls -1
read -p 'What directory? ' -e FOLDER
cd ..

echo "Do you want to delete the images afterwards?"
while true;
do
    read -p "(yes/no): " DEL
    if [ "$DEL" = "yes" ] || [ "$DEL" = "no" ];
    then
        break
    else
        echo "Invalid input."
    fi
done

cp -f input/$FOLDER*.L2D input/$FOLDER*.L3D ./

for FILE in $(find input/$FOLDER*.ini)
do
    echo "Processing $FILE"
    ./engine $FILE
done

eog input/$FOLDER*.bmp input/$FOLDER*.png

rm -f *.L3D *.L2D

if [ "$DEL" == "yes" ];
then
    rm input/$FOLDER*.bmp
fi
