      module fgraphics_mod

      interface
         subroutine C_set_graphics_mysize(ncells) &
              bind(C,NAME="set_graphics_mysize")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value :: ncells
         end subroutine C_set_graphics_mysize

         subroutine C_set_graphics_viewmode(view_mode) &
              bind(C,NAME="set_graphics_viewmode")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value :: view_mode
         end subroutine C_set_graphics_viewmode
         
         subroutine C_set_graphics_cell_coordinates(x, dx, y, dy) &
              bind(C,NAME="set_graphics_cell_coordinates_float")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value :: x, dx, y, dy
         end subroutine C_set_graphics_cell_coordinates
         
         subroutine C_set_graphics_cell_data(data) &
              bind(C,NAME="set_graphics_cell_data_float")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value :: data
         end subroutine C_set_graphics_cell_data
        
         subroutine C_set_graphics_outline(outline) &
              bind(C,NAME="set_graphics_outline")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value :: outline
         end subroutine C_set_graphics_outline
!!!set_graphics_cellproc(&mesh->proc[0]) - need subroutine
!write graphics
!next graphics


         subroutine C_set_graphics_window(xmin, xmax, ymin, ymax) &
              bind(C,NAME="set_graphics_window")
            use iso_c_binding, only : C_FLOAT
            implicit none
            real(C_FLOAT), value :: xmin, xmax, ymin, ymax
         end subroutine C_set_graphics_window

!!! Double check on init and terminate
         subroutine C_terminate_init_graphics_output(argc, argv, title) &
              bind(C,NAME="terminate_graphics_output")
            use iso_c_binding, only : C_PTR, C_CHAR
            implicit none
            type(C_PTR), value :: argc
            type(C_PTR), value :: argv
            character(C_CHAR) :: title
         end subroutine C_terminate_init_graphics_output

         subroutine C_init_graphics_output(argc, argv, title) &
              bind(C,NAME="init_graphics_output")
            use iso_c_binding, only : C_PTR, C_CHAR
            implicit none
            type(C_PTR), value :: argc
            type(C_PTR), value :: argv
            character(C_CHAR) :: title
         end subroutine C_init_graphics_output

         subroutine C_write_graphics_info(graph_num, ncycle, simTime,&
                         rollback_img, rollback_num)&
              bind(C,NAME="write_graphics_info")
            use iso_c_binding, only : C_INT, C_DOUBLE
            implicit none
            integer(C_INT), value :: graph_num, ncycle, rollback_img, rollback_num 
            real(C_DOUBLE), value :: simTIME 
         end subroutine C_write_graphics_info

!         subroutine C_draw_scene() &
!              bind(C,NAME="draw_scene")
!         end subroutine C_draw_scene
      
 

      end interface

      end module fgraphics_mod
