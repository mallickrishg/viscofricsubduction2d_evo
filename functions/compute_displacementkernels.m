function GF_d = compute_displacementkernels(ox,rcv,shz,src)

% rcv
[Gd] = fgeom.computedisplacementGFs(rcv,ox); 
GF_d.rcv.Gh = Gd(1:2:end,:);
GF_d.rcv.Gz = Gd(2:2:end,:);

% shz
[Gd] = fgeom.computedisplacementGFs(shz,ox);
GF_d.shz.Gh = Gd(1:2:end,:);
GF_d.shz.Gz = Gd(2:2:end,:);

% deep source
[Gd] = fgeom.computedisplacementGFs(src,ox);
GF_d.src.Gh = Gd(1:2:end,:);
GF_d.src.Gz = Gd(2:2:end,:);

end