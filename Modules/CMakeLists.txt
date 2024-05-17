art_make( 
		MODULE_LIBRARIES 
                ${ROOT_BASIC_LIB_LIST}
		${ROOT_GEOM}
		${ROOT_XMLIO}
		${ROOT_GDML}
                ${ROOT_GENVECTOR}
		${FHICLCPP}
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${ART_FRAMEWORK_SERVICES_OPTIONAL}
		${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                cetlib cetlib_except
		nusimdata_SimulationBase
		lardataobj_RecoBase
		lardata_RecoObjects
		art_Persistency_Common canvas
		art_Persistency_Provenance canvas
		art_Utilities canvas
)
install_source()
