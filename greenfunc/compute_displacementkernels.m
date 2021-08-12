function GF_d = compute_displacementkernels(ox,rcv,shz,src)

% rcv
[~,Gd] = rcv.displacementKernels(ox,3);
GF_d.rcv.Gh = Gd(1:3:end,:);
GF_d.rcv.Gz = Gd(3:3:end,:);

% shz
[~,Gd] = shz.displacementKernels(ox,3);
GF_d.shz.Gh = Gd(1:3:end,:);
GF_d.shz.Gz = Gd(3:3:end,:);

% deep source
[~,Gd] = src.displacementKernels(ox,3);
GF_d.src.Gh = Gd(1:3:end,:);
GF_d.src.Gz = Gd(3:3:end,:);

end