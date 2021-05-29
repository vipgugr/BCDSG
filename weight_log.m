function w = weight_log(nu)
    val=eps+(abs(nu)+eps).*abs(nu);
    w=1./val;
end