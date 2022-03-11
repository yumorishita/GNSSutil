## Path
export POSDIR="$HOME/GNSS" # May change 

export GNSSutil_PATH="$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)"
export PATH="$GNSSutil_PATH/bin:$PATH"
export PYTHONPATH="$GNSSutil_PATH/GNSSutil_lib:$PYTHONPATH"

