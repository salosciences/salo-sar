if test -d "Retrieved Data";then
    echo "Directory exists"
else
    mkdir "Retrieved Data"
fi

cd "Retrieved Data"

if test -d "wa_v2"; then
    echo "Directory exists"
else
    wget https://gis1.servirglobal.net/TrainingMaterials/SAR/Data_ch3_wa.zip
    unzip Data_ch3_wa.zip
fi
