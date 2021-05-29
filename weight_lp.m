function w = weight_lp(nu,p)
    val=eps+abs(nu).^(2-p);
    w=1./val;
end