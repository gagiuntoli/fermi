function [value] = grad_shape_3_4(i,d,x,y,z)

    switch (i)
      case 1
      % 1-x-y-z
        switch (d)
           case 1
               value = -1 ;
               return
           case 2
              value =  -1 ;
              return
           case 3
              value =  -1 ;
              return
        endswitch
      case 2
      % x
        switch (d)
           case 1
               value = 1 ;
               return
           case 2
              value =  0 ;
              return
           case 3
              value =  0 ;
              return
        endswitch
      case 3
      % y
        switch (d)
           case 1
               value = 0 ;
               return
           case 2
              value =  1 ;
              return
           case 3
              value =  0 ;
              return
        endswitch
      case 4
      % z
        switch (d)
           case 1
               value = 0 ;
               return
           case 2
              value =  0 ;
              return
           case 3
              value =  1 ;
              return
        endswitch       
        
    endswitch

endfunction