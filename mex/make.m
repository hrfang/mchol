% This make.m mimics matlab/make.m in LIBSVM and LIBLINEAR for compiling mex files.
% It has been tested in OCTAVE 3.8.1 under Ubuntu Linux.
% It possibly also works for MATLAB, OCTAVE under Unix, Windows, Mac.

try
    Type = ver;
    % this part is for OCTAVE
    if (strcmp(Type(1).Name, 'Octave') == 1)
        mex gmw81.cpp mldl.cpp ../source/choltool.cpp ../source/gmw.cpp
        mex gmw1.cpp  mldl.cpp ../source/choltool.cpp ../source/gmw.cpp
        mex gmw2.cpp  mldl.cpp ../source/choltool.cpp ../source/gmw.cpp
        mex se90.cpp  mldl.cpp ../source/choltool.cpp ../source/se.cpp
        mex se99.cpp  mldl.cpp ../source/choltool.cpp ../source/se.cpp
        mex se1.cpp   mldl.cpp ../source/choltool.cpp ../source/se.cpp
    % this part is for MATLAB
    % add "-largeArrayDims" on 64-bit machines of MATLAB
    else
        mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims gmw81.cpp mldl.cpp ../source/choltool.cpp ../source/gmw.cpp
        mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims gmw1.cpp  mldl.cpp ../source/choltool.cpp ../source/gmw.cpp
        mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims gmw2.cpp  mldl.cpp ../source/choltool.cpp ../source/gmw.cpp
        mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims se90.cpp  mldl.cpp ../source/choltool.cpp ../source/se.cpp
        mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims se99.cpp  mldl.cpp ../source/choltool.cpp ../source/se.cpp
        mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims se1.cpp   mldl.cpp ../source/choltool.cpp ../source/se.cpp
    end
catch
    fprintf('Oops.. make.m fails..\n');
end
