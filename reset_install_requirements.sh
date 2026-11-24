#!/bin/zsh
# Build the Basilisk source committed alongside this project.
set -eu

if [[ $# == 1 && "$1" == "--hard" ]]; then
    hard=true
elif [[ $# == 0 ]]; then
    hard=false
else
    printf 'Usage: %s [--hard]\n' "$0" >&2
    exit 2
fi

root=$(cd -- "$(dirname -- "$0")" && pwd -P)
basilisk="$root/basilisk/src"
# The committed tree carries local fixes, such as the bound on field names in
# restore(). A tree fetched from upstream by an older copy of this script lacks
# them, so refuse to build it.
if [[ ! -f "$basilisk/Makefile" || ! -f "$basilisk/qcc.c" ]] ||
   ! grep -Fq 'A dump is external data' "$basilisk/output.h" 2>/dev/null; then
    printf '%s\n' 'Restore the vendored basilisk/ directory from this project checkout before building.' >&2
    exit 1
fi

(
    cd -- "$basilisk"
    if [[ "$OSTYPE" == darwin* ]]; then
        platform_config=config.osx
    else
        platform_config=config.gcc
    fi
    # The checkout carries the config link and binaries of the machine that
    # committed them. Keep a custom config; point a stock one at this platform.
    current_config=$(readlink config 2>/dev/null || true)
    current_config=${current_config#./}
    if [[ ! -e config && ! -L config ]] ||
       [[ "$current_config" == config.osx || "$current_config" == config.gcc ]]; then
        ln -sfn "$platform_config" config
    fi
    # Rebuild from source for --hard, or when the binaries cannot run here, or
    # were built somewhere else: qcc looks for its headers in the directory it
    # was built in, which for the committed binaries is another machine's.
    # "make -B" rebuilds ast/interpreter from the wrong directory and fails.
    if [[ "$hard" == true ]] || ! ./qcc --version >/dev/null 2>&1 ||
       ! grep -Fq -- "$basilisk" qcc; then
        make clean
        rm -f ast/grammar  # not removed by "make clean"; make rebuilds it
    fi
    # Serial: the Makefile does not order qcc after the ast sub-build.
    make -j1
)

# qcc makes a temporary directory where it runs, so check it from the source tree.
(cd -- "$basilisk" && ./qcc --version)
# Single quotes keep the path literal for any shell that reads this file,
# including the /bin/sh the simulation Makefiles use.
quote="'"
{
    printf "export BASILISK='%s'\n" "${basilisk//$quote/$quote\\$quote$quote}"
    printf '%s\n' 'export PATH="$BASILISK:$PATH"'
} > "$root/.project_config"
