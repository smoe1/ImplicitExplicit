function q = qexact_slices(time,x,y,z,mdir,params)

    switch(mdir)
      case 1
        [m1,m2]=size(x);        
        q = zeros(m1,m2);
        for i=1:m1
            for j=1:m2
                q(i,j) = qexfunc(time,x(i,j),y(i,j),z,params);
            end
        end
      case 2
        [m1,m2]=size(x);
        q = zeros(m1,m2);
        for i=1:m1
            for j=1:m2
                q(i,j) = qexfunc(time,x(i,j),y,z(i,j),params);
            end
        end
      case 3
        [m1,m2]=size(y);
        q = zeros(m1,m2);
        for i=1:m1
            for j=1:m2
                q(i,j) = qexfunc(time,x,y(i,j),z(i,j),params);        
            end
        end
    end

end