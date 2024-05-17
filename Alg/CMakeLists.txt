art_make( 

      LIBRARY_NAME     uboonecode_uboone_FlexiPID_Alg

      LIB_LIBRARIES
      ${ROOT_BASIC_LIB_LIST}
      ${FHICLCPP}
      cetlib cetlib_except
      art_Persistency_Common canvas
      art_Persistency_Provenance canvas
      art_Utilities canvas
      ${ART_FRAMEWORK_CORE}
      ${ART_FRAMEWORK_PRINCIPAL}
      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
)

install_headers()
install_source()
