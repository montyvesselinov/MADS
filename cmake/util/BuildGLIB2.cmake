build_glib2()
{
    version="$1"
    maj_min_version="${version%.*}" #Drop micro#

    if [ -e $DEPLOYDIR/lib/glib-2.0 ]; then
    echo "glib2 already installed. not building"
    return
    fi

    echo "Building glib2 $version..."
    cd "$BASEDIR"/src
    rm -rf "glib-$version"
    if [ ! -f "glib-$version.tar.xz" ]; then
    curl --insecure -LO "http://ftp.gnome.org/pub/gnome/sources/glib/$maj_min_version/glib-$version.tar.xz"
    fi
    tar xJf "glib-$version.tar.xz"
    cd "glib-$version"
    ./configure --disable-gtk-doc --disable-man --prefix="$DEPLOYDIR" CFLAGS="-I$DEPLOYDIR/include" LDFLAGS="-L$DEPLOYDIR/lib"
    make -j$NUMCPU
    make install
}
