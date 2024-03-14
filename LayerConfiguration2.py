import sys
import re

inputFile = open(sys.argv[1], 'r')
numberset_y = []
numberset_x = []
numberset_z = []
nmodules = ""
scintThickness = ""
scintCasing = ""
wall_gap = ""
wall_gap2 = ""
airGap = ""
layerSpacing = ""
while True:
    line = inputFile.readline()
    if "//" in line:
        continue
    if "if (new_position.y() >=" in line: # This is the place parameters are specified 
        numbers_y = re.findall(r"[-+]?(?:\d*\.*\d+)", line) #regex expression extracts decimal numbers
        numberset_y.append(numbers_y)
    elif "if (new_position.x() >=" in line:
        numbers_x = re.findall(r"[-+]?(?:\d*\.*\d+)", line)
        numberset_x.append(numbers_x)
    elif "if (new_position.z() >=" in line:
        numbers_z = re.findall(r"[-+]?(?:\d*\.*\d+)", line)
        numberset_z.append(numbers_z)
    elif "constexpr auto wall_gap " in line:
        wall_gap = float(re.findall(r"[-+]?(?:\d*\.*\d+)", line)[0]) # in meters
    elif "constexpr auto wall_gap2" in line:
        wall_gap2 = float(re.findall(r"[-+]?(?:\d*\.*\d+)", line)[1]) # avoiding digit in variable name
    elif "constexpr int NMODULES" in line:
        nmodules = re.findall(r"[-+]?(?:\d*\.*\d+)", line)[0]
    elif "constexpr auto scintillator_height =" in line:
        scintThickness = float(re.findall(r"[-+]?(?:\d*\.*\d+)", line)[0]) # in meters
    elif "constexpr auto scintillator_casing_thickness = " in line:
        scintCasing = float(re.findall(r"[-+]?(?:\d*\.*\d+)", line)[0]) # in meters
    elif "constexpr auto air_gap = " in line:
        airGap = float(re.findall(r"[-+]?(?:\d*\.*\d+)", line)[0]) # in meters 
    elif "constexpr auto layer_spacing =" in line: 
        layerSpacing = float(re.findall(r"[-+]?(?:\d*\.*\d+)", line)[0]) # in meters 
    elif "if (_rotation == 1)" in line or line == "":
        break
inputFile.close()

outputFile = open(sys.argv[2], 'r')
lines = outputFile.readlines()
outputFile.close()
outputFile = open(sys.argv[2], 'w')
linenum = 0
while linenum < len(lines):
    if "const std::vector<std::vector<double>> LAYERS_Y=" in lines[linenum]: # Check if this is the vector to change
        outputFile.write("const std::vector<std::vector<double>> LAYERS_Y={{" + numberset_y[0][0] + '*cm,' + numberset_y[0][1] + "*cm},\n")
        for i in range(1, len(numberset_y)):
            outputFile.write('{' + numberset_y[i][0] + '*cm,' + numberset_y[i][1] + "*cm},\n")
        while "*cm" in lines[linenum + 1]: # Skip over previous values, don't write them.
            linenum += 1
        if "};" not in lines[linenum + 1]:
            outputFile.write("};\n") # closing the vector of vectors 
    elif "const std::vector<std::vector<double>> MODULE_X = {" in lines[linenum]:
        outputFile.write("const std::vector<std::vector<double>> MODULE_X = {{" + numberset_x[0][0] + '*cm,' + numberset_x[0][1] + "*cm},\n")
        for i in range(1, len(numberset_x)):
            outputFile.write('{' + numberset_x[i][0] + '*cm,' + numberset_x[i][1] + "*cm},\n")
        while "*cm" in lines[linenum + 1]:
            linenum += 1
        if "};" not in lines[linenum + 1]:
            outputFile.write("};\n") # closing the vector of vectors 
    elif "const std::vector<std::vector<double>> MODULE_Z = {" in lines[linenum]:
        outputFile.write("const std::vector<std::vector<double>> MODULE_Z = {{" + numberset_z[0][0] + '*cm,' + numberset_z[0][1] + "*cm},\n")
        for i in range(1, len(numberset_z)):
            outputFile.write('{' + numberset_z[i][0] + '*cm,' + numberset_z[i][1] + "*cm},\n")
        while "*cm" in lines[linenum + 1]:
            linenum += 1
        if "};" not in lines[linenum + 1]:
            outputFile.write("};\n") # closing the vector of vectors 
    elif "const int n_modules =" in lines[linenum]:
        outputFile.write("const int n_modules = " + nmodules + ";\n")
    elif "const double scintillator_height = " in lines[linenum]:
        outputFile.write("const double scintillator_height = " + str(scintThickness*100 + 2*scintCasing*100) + "*units::cm;\n")
    elif "const double scintillator_thickness = " in lines[linenum]:
        outputFile.write("const double scintillator_thickness = " + str(scintThickness*100) + "*units::cm; //Just the sensitive part\n")
    elif "const double wall_gap = " in lines[linenum]:
        outputFile.write("const double wall_gap = " + str(wall_gap*100) + "*units::cm; //gap on each side of wall\n")
    elif "const double wall_gap2 = " in lines[linenum]:
        outputFile.write("const double wall_gap2 = " + str(wall_gap2*100) + "*units::cm; //gap between two walls\n")
    elif "const double wall_height =" in lines[linenum]:
        outputFile.write("const double wall_height = " + str((airGap + 2*(scintThickness + 2*scintCasing) + layerSpacing)*100) + "*units::cm;\n")
    else:
        outputFile.write(lines[linenum])
    linenum += 1
