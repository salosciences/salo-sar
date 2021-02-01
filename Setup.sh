if test -d "/Retrieved Data/SAR_CH3/";then
    echo "Directory exists"
else
    mkdir "Retrieved Data"
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

if test -d "SAR_CH4";then
    echo "Directory exists"
else
    mkdir "SAR_CH4"
fi

cd "SAR_CH4"

if test -d "test_example_ROIPAC"; then
    echo "Directory exists"
else

fi
