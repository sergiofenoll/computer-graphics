#! /bin/bash

echo 'Available choices:'
cd input/
ls -1
read -p 'What directory? ' -e FOLDER
cd ..

echo "Do you want to open the images afterwards?"
while true;
do
    read -p "(yes/no): " OPEN
    if [ "$OPEN" = "yes" ] || [ "$OPEN" = "no" ];
    then
        break
    else
        echo "Invalid input."
    fi
done

cp input/$FOLDER*.L2D input/$FOLDER*.L3D ./ 2> /dev/null

for FILE in $(find input/$FOLDER*.ini)
do
    echo "Processing $FILE"
    ./engine $FILE
done

rm ./*.L3D ./*.L2D 2> /dev/null

if [ "$OPEN" == "yes" ];
then
    eog input/$FOLDER*.bmp input/$FOLDER*.png
fi
