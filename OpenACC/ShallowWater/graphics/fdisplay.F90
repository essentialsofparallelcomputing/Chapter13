      module fdisplay_mod

      interface
         subroutine C_set_display_mysize(ncells) &
              bind(C,NAME="set_display_mysize")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value :: ncells
         end subroutine C_set_display_mysize

         subroutine C_set_display_window(xmin, xmax, ymin, ymax) &
              bind(C,NAME="set_display_window")
            use iso_c_binding, only : C_FLOAT
            implicit none
            real(C_FLOAT), value :: xmin, xmax, ymin, ymax
         end subroutine C_set_display_window

         subroutine C_set_display_outline(outline) &
              bind(C,NAME="set_display_outline")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value :: outline
         end subroutine C_set_display_outline

         subroutine C_set_display_cell_coordinates(x, dx, y, dy) &
              bind(C,NAME="set_display_cell_coordinates_float")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value :: x, dx, y, dy
         end subroutine C_set_display_cell_coordinates

         !For particle routine
        ! subroutine C_set_number_of_particles_to_display(npart) &
        !   bind(C,NAME="set_number_of_particles_to_display")
        !   use iso_c_binding, only : C_INT
        !   implicit none
        !   integer(C_INT), value :: npart
        !end subroutine C_set_number_of_particles_to_display
  
         subroutine C_set_display_particle_coordinates(xpart,ypart,npart) &
             bind(C,NAME="set_display_particle_coordinates_float")
             use iso_C_binding, only : C_PTR, C_INT
             implicit none
             type(C_PTR), value :: xpart,ypart
             integer(C_INT), value :: npart
         end subroutine C_set_display_particle_coordinates

         subroutine C_set_display_cell_data(data) &
              bind(C,NAME="set_display_cell_data_float")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value :: data
         end subroutine C_set_display_cell_data

         subroutine C_set_display_viewmode(view_mode) &
              bind(C,NAME="set_display_viewmode")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value :: view_mode
         end subroutine C_set_display_viewmode

         subroutine C_init_display(argc, argv, title) &
              bind(C,NAME="init_display")
            use iso_c_binding, only : C_PTR, C_CHAR
            implicit none
            type(C_PTR), value :: argc
            type(C_PTR), value :: argv
            character(C_CHAR) :: title
         end subroutine C_init_display

         subroutine C_draw_scene() &
              bind(C,NAME="draw_scene")
         end subroutine C_draw_scene

         subroutine C_provide_sim_progress(simTime, ncycle) &
              bind(C,NAME="provide_sim_progress")
            use iso_c_binding, only : C_INT, C_DOUBLE
            implicit none
            integer(C_INT), value :: ncycle
            real(C_DOUBLE), value :: simTime
         end subroutine C_provide_sim_progress

      end interface

      end module fdisplay_mod
