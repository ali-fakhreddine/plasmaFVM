#!/bin/bash

# Set the source and module files
MODULE_FILES=(
    "m_precision.f90"
    "m_varspace.f90"
    "m_constants.f90"
    "m_grid.f90"
    "m_dataio.f90"
    "m_init.f90"
    "m_solver.f90"
)
MAIN_FILE="main.f90"
OUTPUT_EXEC="plasmafvm"

# Debug flags
DEBUG_FLAGS="-g -Wall -fbacktrace -fcheck=all -O0"

# Clean up any old .mod or .o files
echo "Cleaning up old files..."
rm -f *.mod *.o 

# Compile each module file
echo "Compiling module files..."
for module_file in "${MODULE_FILES[@]}"; do
    gfortran -c "$module_file"
    if [ $? -ne 0 ]; then
        echo "Error compiling $module_file"
        exit 1
    fi
done

# Compile the main program and link with the module object files
echo "Compiling and linking $MAIN_FILE..."
gfortran $DEBUG_FLAGS -o "$OUTPUT_EXEC" "$MAIN_FILE" *.o

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful! Running the program..."
    # Run the program
    ./$OUTPUT_EXEC
else
    echo "Compilation failed!"
fi
