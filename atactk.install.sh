#the current directory should contain make_cut_matrix.patch, and metrics.py.patch
cur=`pwd`

mkdir git
cd git
git clone https://github.com/ParkerLab/atactk
cd atactk
git checkout 6cd7de0
cd ..
pip install --user ./atactk
cd ..

cd ~/.local/bin
cp $cur/make_cut_matrix.patch .
patch make_cut_matrix < make_cut_matrix.patch
cd ~/.local/lib/python2.7/site-packages/atactk
cp $cur/metrics.py.patch .
patch metrics.py < metrics.py.patch

cd $cur
