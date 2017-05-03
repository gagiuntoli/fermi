function [value] = shape_3_4(i,x,y,z)
      
      switch (i)
          case 1
             % 1-x-y-z
             value = 1-x-y-z;
             return
          case 2
             % x
             value = x;
             return
          case 3
             % y
             value = y;
             return
          case 4
             % z
             value = z;
             return
      endswitch
      return
endfunction