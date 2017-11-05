file_id = H5F.open('data2.h5')
%h5disp('data2.h5')
energy = hdf5read('data2.h5','FES/chi_005/energy');

