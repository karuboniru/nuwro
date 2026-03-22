# cmake/generate_params.cmake
# Called by add_custom_command to build params_all.h from params.xml + params.sed.
# Variables expected (passed via -D):
#   OUTPUT   - path to write params_all.h
#   SED_FILE - path to params.sed
#   XML_FILE - path to params.xml

execute_process(
  COMMAND sed -f "${SED_FILE}" "${XML_FILE}"
  OUTPUT_VARIABLE SED_OUTPUT
  RESULT_VARIABLE SED_RESULT
)

if(NOT SED_RESULT EQUAL 0)
  message(FATAL_ERROR "sed failed (exit ${SED_RESULT}) while generating ${OUTPUT}")
endif()

# Write the macro header followed by the sed-generated PARAM() body.
# Each line produced by sed already ends with a backslash continuation.
file(WRITE "${OUTPUT}" "#define PARAMS_ALL()\\\n${SED_OUTPUT}\n")
