# ISCE Processing for ALOS data...

import os
import numpy as np
import sys


# For data download to succeed, user must have .netrc file in their home directory to store their username and password.
def configure_inputs(outDir, date1, date2, file1, file2):
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

        <property name="do rubbersheetingAzimuth">True</property>
        <property name="do rubbersheetingRange">False</property>
        <property name="do denseoffsets">True</property>
        <property name="do split spectrum">True</property>
        <property name="unwrapper name">snaphu</property> 
        <property name="do dispersive">True</property>
        <property name="dispersive filter kernel x-size">800</property>
        <property name="dispersive filter kernel y-size">800</property>
        <property name="dispersive filter kernel sigma_x">100</property>
        <property name="dispersive filter kernel sigma_y">100</property>
        <property name="dispersive filter kernel rotation">0</property>
        <property name="dispersive filter number of iterations">5</property>
        <property name="dispersive filter mask type">connected_components</property>
        <property name="dispersive filter coherence threshold">0.6</property>


    </component>
    </stripmapApp>"""

    print("writing stripmapApp.xml")
    with open(os.path.join(outDir, "stripmapApp.xml"), "w") as fid:
        fid.write(cmd_stripmap_config)


def ALOS_XML(link1, link2):
    ASF_USER = ""
    ASF_PWD = ""

    home_dir = os.path.join(os.getenv("HOME"), "work")
    process_dir = os.path.join(home_dir, "Tahoe_ALOS1")
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
        else:
            print(zip1, " already exists")
        cmd2 = ("wget  %s --user {0} --password {1}" % link2).format(ASF_USER, ASF_PWD)
        zip2 = os.path.basename(os.path.normpath(link2))  # isolating initial zip file
        if not os.path.exists(os.path.join(data_dir, zip2)):
            os.system(cmd2)
        else:
            print(zip2, " already exists")

    os.system("unzip  %s" % zip1)
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

    configure_inputs(process_dir, date1, date2, zip1[:-4], zip2[:-4])  # Writing reference secondary and stripmap XML
    os.getcwd()

    # os.system("stripmapApp.py stripmapApp.xml --start=startup --end=endup")
