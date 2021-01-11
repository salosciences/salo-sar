if test -d "Retrieved Data/SAR_CH3";then
    echo "Directory exists"
else
    mkdir "Retrieved Data/SAR_CH3"
fi

cd "Retrieved Data/SAR_CH3"

if test -d "wa_v2"; then
    echo "Directory exists"
else
    wget https://gis1.servirglobal.net/TrainingMaterials/SAR/Data_ch3_wa.zip
    unzip Data_ch3_wa.zip
fi

cd ".."

if test -d "SAR_CH2";then
    echo "Directory exists"
else
    mkdir "SAR_CH2"
fi

cd "SAR_CH2"

if test -d "Data_ch2_E2"; then
    echo "Directory exists"
else
    wget https://gis1.servirglobal.net/TrainingMaterials/SAR/Data_ch2_E2.zip
    unzip Data_ch2_E2.zip
fi