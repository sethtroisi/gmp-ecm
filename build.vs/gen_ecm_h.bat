@echo off
echo creating ecm.h from ecm.h.in
echo /* generated from ecm-h.in by gen_ecm_h.bat */>tmp.h

for /f "tokens=1,2*" %%a in (..\ecm.h.in) do (
  if "%%a" EQU "#undef" (
    if "%%b" EQU "ECM_VERSION" (
      echo #define ECM_VERSION "7.0.5-dev">>tmp.h
    )
  ) else echo %%a %%b %%c>>tmp.h
)

call out_copy_rename tmp.h ..\ ecm.h
