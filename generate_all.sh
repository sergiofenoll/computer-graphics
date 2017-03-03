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

cp input/$FOLDER/*.ini input/$FOLDER/*.L2D ./ &>/dev/null
for FILE in $(find *.ini)
do
    echo "Processing $FILE"
    ./engine $FILE
done

eog *.bmp

rm *.ini *.L2D &>/dev/null

if [ "$DEL" == "yes" ];
then
    rm *.bmp
fi
