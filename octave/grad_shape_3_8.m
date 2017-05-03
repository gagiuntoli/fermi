function [value] = grad_shape_3_8(i,d,x,y,z)

    switch (i)
      case 1
      % (1-x)*(1-y)*(1-z)
        switch (d)
           case 1
               value = -1*(1-y)*(1-z) ;
               return
           case 2
              value =  -1*(1-x)*(1-z) ;
              return
           case 3
              value =  -1*(1-x)*(1-y) ;
              return
        endswitch
      case 2
      % x*(1-y)*(1-z)
        switch (d)
           case 1
               value = (1-y)*(1-z) ;
               return
           case 2
              value =  -1*x*(1-z) ;
              return
           case 3
              value =  -1*x*(1-y) ;
              return
        endswitch
      case 3
      % x*y*(1-z)
        switch (d)
           case 1
               value = y*(1-z) ;
               return
           case 2
              value =  x*(1-z) ;
              return
           case 3
              value =  -1*x*y ;
              return
        endswitch
      case 4
      % (1-x)*y*(1-z)
        switch (d)
           case 1
               value = -1*y*(1-z) ;
               return
           case 2
              value =  (1-x)*(1-z) ;
              return
           case 3
              value =  -1*(1-x)*y ;
              return
        endswitch
      case 5
      % (1-x)*(1-y)*z
        switch (d)
           case 1
               value = -1*(1-y)*z ;
               return
           case 2
              value =  -1*(1-x)*z ;
              return
           case 3
              value =  (1-x)*(1-y) ;
              return
        endswitch
      case 6
      % x*(1-y)*z
        switch (d)
           case 1
               value = (1-y)*z ;
               return
           case 2
              value =  -1*x*z ;
              return
           case 3
              value =  x*(1-y) ;
              return
        endswitch
      case 7
      % x*y*z
        switch (d)
           case 1
               value = y*z ;
               return
           case 2
              value =  x*z ;
              return
           case 3
              value =  x*y ;
              return
        endswitch
      case 8
      % (1-x)*y*z
        switch (d)
           case 1
               value = -1*y*z ;
               return
           case 2
              value =  (1-x)*z ;
              return
           case 3
              value =  (1-x)*y ;
              return
        endswitch      
        
        
    endswitch

endfunction