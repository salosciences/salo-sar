# ISCE Processing for ALOS data...

import os
import numpy as np
import sys


# For data download to succeed, user must have .netrc file in their home directory to store their username and password.
def Write_isce_xml(outDir):
    isce = '''<?xml version="1.0" encoding="UTF-8"?>

    <!-- NOTE: tag/attribute names must be in lower case -->
    <!--
    This file can be used to set the steps to use in the flow.
    The default values for these is set to False in the code.
    You can use this file as defaults in a directory pointed to
    by an environment variable named $ISCEDB as your defaults,
    or you can put it in the processing directory with more
    the defaults for the current working directory.  If you have
    both, then the one in the $ISCEDB directory is loaded first
    and then the one in the local working directory is loaded.
    Additional configuration parameters can be given in another
    file named isceApp.xml that will be loaded automatically
    without putting it on the command line or else you can give
    it any desired name and put it on the command line.  Exmaples
    in this directory named isceappUAVSAR_Stack.xml and
    isceappALOS_pol.xml give examples."
    -->

    <isceApp>
    <component name="isce">
        <!-- Processors to run: True/False -->
        <property name="do preprocess">True</property>
        <property name="do verifyDEM">True</property>
        <property name="do pulsetiming">True</property>
        <property name="do estimateheights">True</property>
        <property name="do mocomppath">True</property>
        <property name="do orbit2sch">True</property>
        <property name="do updatepreprocinfo">True</property>
        <property name="do formslc">True</property>
        <property name="do multilookslc">True</property>
        <property name="do filterslc">False</property>
        <property name="do polarimetric correction">False</property>
        <property name="do calculate FR">False</property>
        <property name="do FR to TEC">False</property>
        <property name="do TEC to phase">False</property>
    <!--    <property name="do geocodeslc">False</property>  -->
    <!--    <property name="do offsetprf">False</property> -->
        <property name="do outliers1">False</property>-->
        <property name="do prepareresamps">True</property>-->
        <property name="do resamp">False</property>-->
        <property name="do resamp image">False</property>-->
        <property name="do crossmul">True</property>
        <property name="do mocomp baseline">True</property>
        <property name="do set topoint1">True</property>
        <property name="do topo">True</property>
        <property name="do shadecpx2rg">True</property>
    <!--    <property name="do rgoffset">True</property> -->
        <property name="do rg outliers2">True</property>
        <property name="do resamp only">True</property>
        <property name="do set topoint2">True</property>
        <property name="do correct">True</property>
        <property name="do coherence">True</property>
        <property name="do filter interferogram">True</property>
    <!--    <property name="unwrap">True</property> -->
        <property name="do unwrap">False</property>
        <property name="do geocode">True</property>
    </component>
    </isceApp>'''

    print("writing Reference.xml")
    with open(os.path.join(outDir, "Isce.xml"), "w") as fid:
        fid.write(isce)

def SENT_configure_inputs(outDir, date1, date2, file1, file2):
    cmd_reference_config = '''<?xml version="1.0" encoding="UTF-8"?>
    <component>
    <property name="safe">
        <value>data/S1A_S1_SLC__1SSV_20170122T234204_20170122T234233_014951_01867B_A254.zip</value>
    </property>
    <property name="orbit directory">
        <value>/home/s1/orbits/poeorb</value>
    </property>
    <property name="OUTPUT">
        <value>20170122</value>
    </property>
    </component>
    <!-- This is for Strip Map Sentinel SLC data only -->
    <!-- Shown above is the bare minimum input XML file for Sentinel1A sensor. The reference and secondary catalog xml files have the same structure.'''
    
    print("writing Reference.xml")
    with open(os.path.join(outDir, "Reference.xml"), "w") as fid:
        fid.write(cmd_reference_config)

    cmd_secondary_config = '''<?xml version="1.0" encoding="UTF-8"?>
    <component>
    <property name="safe">
        <value>data/S1A_S1_SLC__1SSV_20170122T234204_20170122T234233_014951_01867B_A254.zip</value>
    </property>
    <property name="orbit directory">
        <value>/home/s1/orbits/poeorb</value>
    </property>
    <property name="OUTPUT">
        <value>20170122</value>
    </property>
    </component>
    <!-- This is for Strip Map Sentinel SLC data only -->
    <!-- Shown above is the bare minimum input XML file for Sentinel1A sensor. The reference and secondary catalog xml files have the same structure.'''
    
    print("writing stripmapApp.xml")
    with open(os.path.join(outDir, "Secondary.xml"), "w") as fid:
        fid.write(cmd_secondary_config)

    cmd_reference_config = '''<?xml version="1.0" encoding="UTF-8"?>
        <insarApp>
        <component name="insar">
            <property  name="Sensor name">ALOS</property>
            <component name="reference">
                <catalog>reference.xml</catalog>
            </component>
            <component name="secondary">
                <catalog>secondary.xml</catalog>
            </component>
        </component>
    </insarApp>'''

    print("writing stripmapApp.xml")
    with open(os.path.join(outDir, "stripmapApp.xml"), "w") as fid:
        fid.write(cmd_reference_config)


def ALOS_configure_inputs(outDir, date1, date2, file1, file2):
    """Write Configuration files for ISCE2 stripmapApp to process NISAR sample products"""
    cmd_reference_config = """<component>
        <property name="IMAGEFILE">
        <value>[data/{0}/{1}/IMG-HH-{2}-H1.0__A]</value>
        </property>
        <property name="LEADERFILE">
            <value>[data/{0}/{1}/LED-{2}-H1.0__A]</value>
        </property>
        <property name="OUTPUT">
            <value>{0}</value>
        </property>
    </component>""".format(
        date1, file1, file1[:-5]
    )

    print("writing reference.xml")
    with open(os.path.join(outDir, "reference.xml"), "w") as fid:
        fid.write(cmd_reference_config)

    cmd_secondary_config = """<component>
        <property name="IMAGEFILE">
            <value>[data/{0}/{1}/IMG-HH-{2}-H1.0__A]</value>
        </property>
        <property name="LEADERFILE">
            <value>[data/{0}/{1}/LED-{2}-H1.0__A]</value>
        </property>
        <property name="OUTPUT">
            <value>{0}</value>
        </property>
    </component>""".format(
        date2, file2, file2[:-5]
    )

    print("writing secondary.xml")
    with open(os.path.join(outDir, "secondary.xml"), "w") as fid:
        fid.write(cmd_secondary_config)

    cmd_stripmap_config = """<?xml version="1.0" encoding="UTF-8"?>
    <stripmapApp>
    <component name="insar">
        <property name="sensor name">ALOS</property>
        <component name="reference">
            <catalog>reference.xml</catalog>
        </component>
        <component name="secondary">
            <catalog>secondary.xml</catalog>
        </component>
    </component>
    </stripmapApp>"""

    print("writing stripmapApp.xml")
    with open(os.path.join(outDir, "stripmapApp.xml"), "w") as fid:
        fid.write(cmd_stripmap_config)


def XML(Sensor, link1, link2):
    ASF_USER = ""
    ASF_PWD = ""

    home_dir = os.path.join(os.getenv("HOME"), "work")
    process_dir = os.path.join(home_dir, "Tahoe_{0}").format(Sensor)
    data_dir = os.path.join(process_dir, "data")

    # If ASF USER/PWD is not provided read it from the .netrc file
    if (len(ASF_PWD) == 0 | len(ASF_USER) == 0) & (os.path.exists(os.path.join(os.getenv("HOME"), ".netrc"))):
        netrc_path = os.path.join(os.getenv("HOME"), ".netrc")
        count = len(open(netrc_path).readlines())
        if count == 1:
            file = open(os.path.join(os.getenv("HOME"), ".netrc"), "r")
            contents = file.read().split(" ")
            ASF_USER = contents[3]
            ASF_PWD = contents[5]
            file.close()
        else:
            ASF_USER = np.loadtxt(os.path.join(os.getenv("HOME"), ".netrc"), skiprows=1, usecols=1, dtype=str)[0]
            ASF_PWD = np.loadtxt(os.path.join(os.getenv("HOME"), ".netrc"), skiprows=1, usecols=1, dtype=str)[1]

    if (len(ASF_PWD) == 0 | len(ASF_USER) == 0) & (not os.path.exists(os.path.join(os.getenv("HOME"), ".netrc"))):
        print("WARNING: The ASF User and Password needs to be included in ~/.netrc file.")
        print(" The ~/.netrc file does not exist or has not been setup properly")
        sys.exit("System Exiting until ~/.netrc file setup properly")

    # Now lets check if processing and data directories exist
    # If they do not we will make the directories

    if not os.path.exists(process_dir):
        print("create ", process_dir)
        os.makedirs(process_dir)
    else:
        print(process_dir, " exists")

    if not os.path.exists(data_dir):
        print("create ", data_dir)
        os.makedirs(data_dir)
    else:
        print(data_dir, " exists")

    os.chdir(data_dir)
    # Now lets download the data to the data directory
    if len(link1) == 0 | len(link2) == 0:
        sys.exit("Data links do not exist")
    else:
        cmd1 = ("wget  %s --user {0} --password {1}" % link1).format(ASF_USER, ASF_PWD)
        zip1 = os.path.basename(os.path.normpath(link1))  # isolating initial zip file
        if not os.path.exists(os.path.join(data_dir, zip1)):
            os.system(cmd1)
            os.system("unzip  %s" % zip1)
        else:
            print(zip1, " already exists")
        cmd2 = ("wget  %s --user {0} --password {1}" % link2).format(ASF_USER, ASF_PWD)
        zip2 = os.path.basename(os.path.normpath(link2))  # isolating initial zip file
        if not os.path.exists(os.path.join(data_dir, zip2)):
            os.system(cmd2)
        else:
            print(zip2, " already exists")
            os.system("unzip  %s" % zip2)
    
    SCDT1 = "{0}/{1}".format(zip1[:-4], "%s.l0.workreport" % zip1[:-9])
    with open(SCDT1, "r") as myfile:
        data = myfile.readlines()
        date1 = data[6]
        date1 = date1[38:46]
    print(date1)

    SCDT2 = "{0}/{1}".format(zip2[:-4], "%s.l0.workreport" % zip2[:-9])
    with open(SCDT2, "r") as myfile:
        data = myfile.readlines()
        date2 = data[6]
        date2 = date2[38:46]
    print(date2)
    if not os.path.exists(os.path.join(data_dir, date1)):
        os.mkdir(date1)
    if not os.path.exists(os.path.join(data_dir, date2)):
        os.mkdir(date2)

    os.system("mv {0} {1}".format(zip1[:-4], date1))
    os.system("mv {0} {1}".format(zip2[:-4], date2))

    os.chdir(process_dir)

    if Sensor == "ALOS1":
        Write_isce_xml(process_dir)
        ALOS_configure_inputs(process_dir, date1, date2, zip1[:-4], zip2[:-4])  # Writing reference secondary and stripmap XML
    elif Sensor == "Sentinel1":
        sys.exit("Sentinel 1 not configured yet")

    # os.system("stripmapApp.py stripmapApp.xml --start=startup --end=endup")
