home_directory=$(pwd)

# python3 prep_files.py
working_directory=$(python3 prep_files.py) # capturing working directory printed from prep_files.py

cd $working_directory # ensures .sim files end up in the right place
source runCosima.sh
find . -type f -name "*.gz" -exec gunzip {} +
source runRevan.sh

cd $home_directory
python3 get_figure_of_merit_plots.py $working_directory
