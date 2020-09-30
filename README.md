# Tinker: Software Tools for Molecular Design

<H2><B>Installation on cluster</B></H2>

First you will need to download and build FFTW3. Once you have downloaded and extracted compile and install it using (replace $USER):
<ol>
  <li>./configure --prefix=/home/$USER</li>
  <li>make</li>
  <li>make install</li>
</ol>

Now add the following line to your .bashrc (or just run it before building, idk):
export PKG_CONFIG_PATH=/home/cuzzocrea/programs/fftw-3.3.8/build/

Now download this repository and do the following:
<ul>
  <li> cd Tinker </li>
  <li> cp cmake/CMakeLists.txt source/CMakeLists.txt </li>
  <li> mkdir build </li>
  <li> cd build </li>
  <li> cmake -DCMAKE_INSTALL_PREFIX=/home/$USER ../source </li>
  <li> make </li>
  <li> make install </li>
</ul>
