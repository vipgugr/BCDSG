function c = getfilters(filtersetname)

    switch filtersetname
        
        case 'none'
            c{1} = [1];
            
        case 'fohv'
            c{1}=[0 1 -1]/sqrt(2.0);
            c{2}=[0 1 -1]'/sqrt(2.0);
            
        case 'fo'
            c{1}=[0 1 -1]/sqrt(2.0);
            c{2}=[0 1 -1]'/sqrt(2.0);
            c{3}=[0 0 0;0 1 0;0 0 -1]/sqrt(2.0);
            c{4}=[0 0 -1;0 1 0;0 0 0]/sqrt(2.0);
            
        otherwise
            c{1} = [1];
            
    end
end