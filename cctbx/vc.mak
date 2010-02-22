CPPFLAGS= -O2 -MT -YX -EHsc -I"C:\Program Files\Microsoft Visual C++ Toolkit 2003\include" -I"C:\Dev\Fox\cctbx\cctbx\include" -I"cctbx\include" -I"scitbx\include" -I"."
#/LIBPATH:"E:\Program Files\Microsoft SDK\Lib\IA64"
LINKFLAGS= /NOLOGO /LIBPATH:"C:\Program Files\Microsoft Visual C++ Toolkit 2003\lib"
LINKLIBS=libc.lib

OBJ_ELTBX= cctbx\eltbx\basic.obj cctbx\eltbx\icsd_radii.obj cctbx\eltbx\tiny_pse.obj cctbx\eltbx\fp_fdp.obj \
           cctbx\eltbx\neutron.obj cctbx\eltbx\wavelengths.obj \
           cctbx\eltbx\xray_scattering\it1992.obj cctbx\eltbx\xray_scattering\n_gaussian.obj \
           cctbx\eltbx\xray_scattering\n_gaussian_raw.obj cctbx\eltbx\xray_scattering\wk1995.obj \
           cctbx\eltbx\henke.obj cctbx\eltbx\henke_tables_25_36.obj  cctbx\eltbx\henke_tables_61_72.obj  \
           cctbx\eltbx\henke_tables_01_12.obj cctbx\eltbx\henke_tables_37_48.obj cctbx\eltbx\henke_tables_73_84.obj  \
           cctbx\eltbx\henke_tables_13_24.obj cctbx\eltbx\henke_tables_49_60.obj cctbx\eltbx\henke_tables_85_92.obj \
           cctbx\eltbx\sasaki.obj cctbx\eltbx\sasaki_tables_25_36.obj cctbx\eltbx\sasaki_tables_61_72.obj \
           cctbx\eltbx\sasaki_tables_01_12.obj cctbx\eltbx\sasaki_tables_37_48.obj cctbx\eltbx\sasaki_tables_73_82.obj \
           cctbx\eltbx\sasaki_tables_13_24.obj cctbx\eltbx\sasaki_tables_49_60.obj

OBJ_SGTBX = cctbx\sgtbx\bricks.obj cctbx\sgtbx\miller.obj cctbx\sgtbx\select_generators.obj \
            cctbx\sgtbx\tr_group.obj cctbx\sgtbx\change_of_basis_op.obj cctbx\sgtbx\reciprocal_space_asu.obj \
            cctbx\sgtbx\seminvariant.obj cctbx\sgtbx\tr_vec.obj \
            cctbx\sgtbx\find_affine.obj cctbx\sgtbx\reciprocal_space_ref_asu.obj cctbx\sgtbx\site_symmetry.obj cctbx\sgtbx\utils.obj \
            cctbx\sgtbx\group_codes.obj cctbx\sgtbx\rot_mx.obj  cctbx\sgtbx\space_group.obj cctbx\sgtbx\wyckoff.obj \
            cctbx\sgtbx\hall_in.obj cctbx\sgtbx\rot_mx_info.obj cctbx\sgtbx\space_group_type.obj \
            cctbx\sgtbx\lattice_symmetry.obj cctbx\sgtbx\row_echelon_solve.obj cctbx\sgtbx\symbols.obj \
            cctbx\sgtbx\lattice_tr.obj cctbx\sgtbx\rt_mx.obj  cctbx\sgtbx\tensor_rank_2.obj \
            cctbx\sgtbx\reference_settings\hall_symbol_table.obj cctbx\sgtbx\reference_settings\normalizer.obj \
            cctbx\sgtbx\reference_settings\matrix_group_code_table.obj cctbx\sgtbx\reference_settings\wyckoff.obj

OBJ_GLOBAL = cctbx\global\error.obj

OBJ_MILLER = cctbx\miller\asu.obj cctbx\miller\index_generator.obj cctbx\miller\match_bijvoet_mates.obj \
             cctbx\miller\sym_equiv.obj cctbx\miller\bins.obj cctbx\miller\index_span.obj cctbx\miller\match_indices.obj

OBJ_UCTBX = cctbx\uctbx\uctbx.obj cctbx\uctbx\spoil_optimization.obj

OBJ = $(OBJ_ELTBX) $(OBJ_SGTBX) $(OBJ_GLOBAL) $(OBJ_MILLER) $(OBJ_UCTBX)

.cpp.obj:
	cl $(CPPFLAGS) -c $*.cpp -Fo$*.obj

libcctbx: $(OBJ)
	lib -nologo -OUT:"libcctbx.lib" $(OBJ)

lib: libcctbx

all: libcctbx

default: lib

tidy:
	del *.obj

clean:
	del $(OBJ) *.lib
