# Abort on errors
set -e

# If on Linux, install dependencies for rust-htslib
case "$(uname -s)" in
    Linux)
        sudo apt install zlib1g-dev libbz2-dev liblzma-dev clang pkg-config
    ;;
    *)
esac