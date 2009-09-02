module forteuchos
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  implicit none   ! Prevent implicit typing

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/epetra/CEpetra_*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

  ! _________________ CommandLineProcessor interface bodies _________________


  ! Original C++ prototype:
  ! ~CommandLineProcessor();
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_Destroy ( CT_Teuchos_CommandLineProcessor_ID_t * selfID );

  subroutine Teuchos_CommandLineProcessor_Destroy ( selfID ) &
        bind(C,name='Teuchos_CommandLineProcessor_Destroy')
    import :: FT_Teuchos_CommandLineProcessor_ID_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t)                                  :: selfID
  end subroutine


  ! Original C++ prototype:
  ! CommandLineProcessor( bool throwExceptions = true ,bool recogniseAllOptions = true ,bool addOutputSetupOptions = false );
  ! CTrilinos prototype:
  ! CT_Teuchos_CommandLineProcessor_ID_t Teuchos_CommandLineProcessor_Create ( boolean throwExceptions, boolean recogniseAllOptions, boolean addOutputSetupOptions );

  type(FT_Teuchos_CommandLineProcessor_ID_t) function Teuchos_CommandLineProcessor_Create ( &
        throwExceptions, recogniseAllOptions, addOutputSetupOptions ) &
        bind(C,name='Teuchos_CommandLineProcessor_Create')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_bool
    
    logical(c_bool)                           ,intent(in)   ,value              :: throwExceptions
    logical(c_bool)                           ,intent(in)   ,value              :: recogniseAllOptions
    logical(c_bool)                           ,intent(in)   ,value              :: addOutputSetupOptions
  end function


  ! Original C++ prototype:
  ! void throwExceptions( const bool & throwExceptions );
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_throwExceptions_set ( CT_Teuchos_CommandLineProcessor_ID_t selfID, boolean throwExceptions );

  subroutine Teuchos_CommandLineProcessor_throwExceptions_set ( selfID, throwExceptions ) &
        bind(C,name='Teuchos_CommandLineProcessor_throwExceptions_set')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_bool
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    logical(c_bool)                           ,intent(in)   ,value              :: throwExceptions
  end subroutine


  ! Original C++ prototype:
  ! bool throwExceptions() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_CommandLineProcessor_throwExceptions_get ( CT_Teuchos_CommandLineProcessor_ID_t selfID );

  logical(c_bool) function Teuchos_CommandLineProcessor_throwExceptions_get ( selfID ) &
        bind(C,name='Teuchos_CommandLineProcessor_throwExceptions_get')
    import :: c_bool ,FT_Teuchos_CommandLineProcessor_ID_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! void recogniseAllOptions( const bool & recogniseAllOptions );
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_recogniseAllOptions_set ( CT_Teuchos_CommandLineProcessor_ID_t selfID, boolean recogniseAllOptions );

  subroutine Teuchos_CommandLineProcessor_recogniseAllOptions_set ( selfID, &
        recogniseAllOptions ) &
        bind(C,name='Teuchos_CommandLineProcessor_recogniseAllOptions_set')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_bool
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    logical(c_bool)                           ,intent(in)   ,value              :: recogniseAllOptions
  end subroutine


  ! Original C++ prototype:
  ! bool recogniseAllOptions() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_CommandLineProcessor_recogniseAllOptions_get ( CT_Teuchos_CommandLineProcessor_ID_t selfID );

  logical(c_bool) function Teuchos_CommandLineProcessor_recogniseAllOptions_get ( selfID ) &
        bind(C,name='Teuchos_CommandLineProcessor_recogniseAllOptions_get')
    import :: c_bool ,FT_Teuchos_CommandLineProcessor_ID_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! void addOutputSetupOptions( const bool &addOutputSetupOptions );
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_addOutputSetupOptions_set ( CT_Teuchos_CommandLineProcessor_ID_t selfID, boolean addOutputSetupOptions );

  subroutine Teuchos_CommandLineProcessor_addOutputSetupOptions_set ( selfID, &
        addOutputSetupOptions ) &
        bind(C,name='Teuchos_CommandLineProcessor_addOutputSetupOptions_set')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_bool
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    logical(c_bool)                           ,intent(in)   ,value              :: addOutputSetupOptions
  end subroutine


  ! Original C++ prototype:
  ! bool addOutputSetupOptions() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_CommandLineProcessor_addOutputSetupOptions_get ( CT_Teuchos_CommandLineProcessor_ID_t selfID );

  logical(c_bool) function Teuchos_CommandLineProcessor_addOutputSetupOptions_get ( selfID ) &
        bind(C,name='Teuchos_CommandLineProcessor_addOutputSetupOptions_get')
    import :: c_bool ,FT_Teuchos_CommandLineProcessor_ID_t
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! void setDocString( const char doc_string[] );
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_setDocString ( CT_Teuchos_CommandLineProcessor_ID_t selfID, const char doc_string[] );

  subroutine Teuchos_CommandLineProcessor_setDocString ( selfID, doc_string ) &
        bind(C,name='Teuchos_CommandLineProcessor_setDocString')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: doc_string
  end subroutine


  ! Original C++ prototype:
  ! void setOption( const char option_true[] ,const char option_false[] ,bool *option_val ,const char documentation[] = NULL );
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_setOption_bool ( CT_Teuchos_CommandLineProcessor_ID_t selfID, const char option_true[], const char option_false[], boolean * option_val, const char documentation[] );

  subroutine Teuchos_CommandLineProcessor_setOption_bool ( selfID, option_true, &
        option_false, option_val, documentation ) &
        bind(C,name='Teuchos_CommandLineProcessor_setOption_bool')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char ,c_bool
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_true
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_false
    logical(c_bool)                                               ,dimension(*) :: option_val
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: documentation
  end subroutine


  ! Original C++ prototype:
  ! void setOption( const char option_name[] ,int *option_val ,const char documentation[] = NULL ,const bool required = false );
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_setOption_int ( CT_Teuchos_CommandLineProcessor_ID_t selfID, const char option_name[], int * option_val, const char documentation[], const boolean required );

  subroutine Teuchos_CommandLineProcessor_setOption_int ( selfID, option_name, option_val, &
        documentation, required ) bind(C,name='Teuchos_CommandLineProcessor_setOption_int')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char ,c_int ,c_bool
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_name
    integer(c_int)                                                ,dimension(*) :: option_val
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: documentation
    logical(c_bool)                           ,intent(in)   ,value              :: required
  end subroutine


  ! Original C++ prototype:
  ! void setOption( const char option_name[] ,double *option_val ,const char documentation[] = NULL ,const bool required = false );
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_setOption_double ( CT_Teuchos_CommandLineProcessor_ID_t selfID, const char option_name[], double * option_val, const char documentation[], const boolean required );

  subroutine Teuchos_CommandLineProcessor_setOption_double ( selfID, option_name, &
        option_val, documentation, required ) &
        bind(C,name='Teuchos_CommandLineProcessor_setOption_double')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char ,c_double ,c_bool
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_name
    real(c_double)                                                ,dimension(*) :: option_val
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: documentation
    logical(c_bool)                           ,intent(in)   ,value              :: required
  end subroutine


  ! Original C++ prototype:
  ! void setOption( const char option_name[] ,std::string *option_val ,const char documentation[] = NULL ,const bool required = false );
  ! CTrilinos prototype:
  ! void Teuchos_CommandLineProcessor_setOption_str ( CT_Teuchos_CommandLineProcessor_ID_t selfID, const char option_name[], char * option_val[], const char documentation[], const boolean required );

  subroutine Teuchos_CommandLineProcessor_setOption_str ( selfID, option_name, option_val, &
        documentation, required ) bind(C,name='Teuchos_CommandLineProcessor_setOption_str')
    import :: FT_Teuchos_CommandLineProcessor_ID_t ,c_char ,c_bool
    
    type(FT_Teuchos_CommandLineProcessor_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: option_name
    character(kind=c_char)                                        ,dimension(*) :: option_val
    character(kind=c_char)                    ,intent(in)         ,dimension(*) :: documentation
    logical(c_bool)                           ,intent(in)   ,value              :: required
  end subroutine


  ! _________________ ParameterList interface bodies _________________


  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Cast ( CTrilinos_Object_ID_t id );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_Cast ( id ) &
        bind(C,name='Teuchos_ParameterList_Cast')
    import :: FT_Teuchos_ParameterList_ID_t ,ForTrilinos_Object_ID_t
    
    type(ForTrilinos_Object_ID_t)      ,intent(in)   ,value              :: id
  end function


  ! CTrilinos prototype:
  ! CTrilinos_Object_ID_t Teuchos_ParameterList_Abstract ( CT_Teuchos_ParameterList_ID_t id );

  type(ForTrilinos_Object_ID_t) function Teuchos_ParameterList_Abstract ( id ) &
        bind(C,name='Teuchos_ParameterList_Abstract')
    import :: ForTrilinos_Object_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: id
  end function


  ! Original C++ prototype:
  ! ParameterList();
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create (  );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_Create (  ) &
        bind(C,name='Teuchos_ParameterList_Create')
    import :: FT_Teuchos_ParameterList_ID_t
    
  end function


  ! Original C++ prototype:
  ! ParameterList(const std::string &name);
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_WithName ( const char name[] );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_Create_WithName ( name ) &
        bind(C,name='Teuchos_ParameterList_Create_WithName')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! ParameterList(const ParameterList& source);
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_Create_FromSource ( CT_Teuchos_ParameterList_ID_t sourceID );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_Create_FromSource ( &
        sourceID ) bind(C,name='Teuchos_ParameterList_Create_FromSource')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: sourceID
  end function


  ! Original C++ prototype:
  ! virtual ~ParameterList();
  ! CTrilinos prototype:
  ! void Teuchos_ParameterList_Destroy ( CT_Teuchos_ParameterList_ID_t * selfID );

  subroutine Teuchos_ParameterList_Destroy ( selfID ) &
        bind(C,name='Teuchos_ParameterList_Destroy')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t)                                  :: selfID
  end subroutine


  ! Original C++ prototype:
  ! ParameterList& setName( const std::string &name );
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setName ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_setName ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_setName')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! ParameterList& operator=(const ParameterList& source);
  ! CTrilinos prototype:
  ! void Teuchos_ParameterList_Assign ( CT_Teuchos_ParameterList_ID_t selfID, CT_Teuchos_ParameterList_ID_t sourceID );

  subroutine Teuchos_ParameterList_Assign ( selfID, sourceID ) &
        bind(C,name='Teuchos_ParameterList_Assign')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: sourceID
  end subroutine


  ! Original C++ prototype:
  ! ParameterList& setParameters(const ParameterList& source);
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParameters ( CT_Teuchos_ParameterList_ID_t selfID, CT_Teuchos_ParameterList_ID_t sourceID );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_setParameters ( selfID, &
        sourceID ) bind(C,name='Teuchos_ParameterList_setParameters')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: sourceID
  end function


  ! Original C++ prototype:
  ! ParameterList& setParametersNotAlreadySet(const ParameterList& source);
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setParametersNotAlreadySet ( CT_Teuchos_ParameterList_ID_t selfID, CT_Teuchos_ParameterList_ID_t sourceID );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_setParametersNotAlreadySet ( &
        selfID, sourceID ) bind(C,name='Teuchos_ParameterList_setParametersNotAlreadySet')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: sourceID
  end function


  ! Original C++ prototype:
  ! ParameterList& disableRecursiveValidation();
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_disableRecursiveValidation ( CT_Teuchos_ParameterList_ID_t selfID );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_disableRecursiveValidation ( &
        selfID ) bind(C,name='Teuchos_ParameterList_disableRecursiveValidation')
    import :: FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! ParameterList& set( std::string const& name, ParameterList const& value, std::string const& docString = "" );
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_set ( CT_Teuchos_ParameterList_ID_t selfID, char const name[], CT_Teuchos_ParameterList_ID_t valueID, char const docString[] );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_set ( selfID, name, &
        valueID, docString ) bind(C,name='Teuchos_ParameterList_set')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: valueID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: docString
  end function


  ! Original C++ prototype:
  ! ParameterList& setEntry(const std::string& name, const ParameterEntry& entry);
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_setEntry ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], CT_Teuchos_ParameterEntry_ID_t entryID );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_setEntry ( selfID, &
        name, entryID ) bind(C,name='Teuchos_ParameterList_setEntry')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: entryID
  end function


  ! Original C++ prototype:
  ! template<typename T> T& get(const std::string& name, T def_value);
  ! CTrilinos prototype:
  ! double Teuchos_ParameterList_get_double_def ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], double def_value );

  real(c_double) function Teuchos_ParameterList_get_double_def ( selfID, name, def_value ) &
        bind(C,name='Teuchos_ParameterList_get_double_def')
    import :: c_double ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    real(c_double)                     ,intent(in)   ,value              :: def_value
  end function


  ! Original C++ prototype:
  ! template<typename T> T& get(const std::string& name, T def_value);
  ! CTrilinos prototype:
  ! int Teuchos_ParameterList_get_int_def ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], int def_value );

  integer(c_int) function Teuchos_ParameterList_get_int_def ( selfID, name, def_value ) &
        bind(C,name='Teuchos_ParameterList_get_int_def')
    import :: c_int ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    integer(c_int)                     ,intent(in)   ,value              :: def_value
  end function


  ! Original C++ prototype:
  ! std::string& get(const std::string& name, char def_value[]);
  ! CTrilinos prototype:
  ! const char * Teuchos_ParameterList_get_char_def ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], char def_value[] );

  character(kind=c_char) function Teuchos_ParameterList_get_char_def ( selfID, name, &
        def_value ) bind(C,name='Teuchos_ParameterList_get_char_def')
    import :: kind=c_char ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    character(kind=c_char)                                 ,dimension(*) :: def_value
  end function


  ! Original C++ prototype:
  ! std::string& get(const std::string& name, const char def_value[]);
  ! CTrilinos prototype:
  ! const char * Teuchos_ParameterList_get_const_char_def ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], const char def_value[] );

  character(kind=c_char) function Teuchos_ParameterList_get_const_char_def ( selfID, name, &
        def_value ) bind(C,name='Teuchos_ParameterList_get_const_char_def')
    import :: kind=c_char ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: def_value
  end function


  ! Original C++ prototype:
  ! template<typename T> T& get(const std::string& name);
  ! CTrilinos prototype:
  ! double Teuchos_ParameterList_get_double ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  real(c_double) function Teuchos_ParameterList_get_double ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_get_double')
    import :: c_double ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> T& get(const std::string& name);
  ! CTrilinos prototype:
  ! int Teuchos_ParameterList_get_int ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  integer(c_int) function Teuchos_ParameterList_get_int ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_get_int')
    import :: c_int ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> const T& get(const std::string& name) const;
  ! CTrilinos prototype:
  ! double Teuchos_ParameterList_get_double_const ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  real(c_double) function Teuchos_ParameterList_get_double_const ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_get_double_const')
    import :: c_double ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> const T& get(const std::string& name) const;
  ! CTrilinos prototype:
  ! int Teuchos_ParameterList_get_int_const ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  integer(c_int) function Teuchos_ParameterList_get_int_const ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_get_int_const')
    import :: c_int ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> inline T* getPtr(const std::string& name);
  ! CTrilinos prototype:
  ! double * Teuchos_ParameterList_getPtr_double ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(c_ptr) function Teuchos_ParameterList_getPtr_double ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_getPtr_double')
    import :: c_ptr ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> inline T* getPtr(const std::string& name);
  ! CTrilinos prototype:
  ! int * Teuchos_ParameterList_getPtr_int ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(c_ptr) function Teuchos_ParameterList_getPtr_int ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_getPtr_int')
    import :: c_ptr ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> inline const T* getPtr(const std::string& name) const;
  ! CTrilinos prototype:
  ! const double * Teuchos_ParameterList_getPtr_double_const ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(c_ptr) function Teuchos_ParameterList_getPtr_double_const ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_getPtr_double_const')
    import :: c_ptr ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> inline const T* getPtr(const std::string& name) const;
  ! CTrilinos prototype:
  ! const int * Teuchos_ParameterList_getPtr_int_const ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(c_ptr) function Teuchos_ParameterList_getPtr_int_const ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_getPtr_int_const')
    import :: c_ptr ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! ParameterEntry& getEntry(const std::string& name);
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(FT_Teuchos_ParameterEntry_ID_t) function Teuchos_ParameterList_getEntry ( selfID, &
        name ) bind(C,name='Teuchos_ParameterList_getEntry')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! inline const ParameterEntry& getEntry(const std::string& name) const;
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntry_const ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(FT_Teuchos_ParameterEntry_ID_t) function Teuchos_ParameterList_getEntry_const ( &
        selfID, name ) bind(C,name='Teuchos_ParameterList_getEntry_const')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! inline ParameterEntry* getEntryPtr(const std::string& name);
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(FT_Teuchos_ParameterEntry_ID_t) function Teuchos_ParameterList_getEntryPtr ( selfID, &
        name ) bind(C,name='Teuchos_ParameterList_getEntryPtr')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! inline const ParameterEntry* getEntryPtr(const std::string& name) const;
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterList_getEntryPtr_const ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(FT_Teuchos_ParameterEntry_ID_t) function Teuchos_ParameterList_getEntryPtr_const ( &
        selfID, name ) bind(C,name='Teuchos_ParameterList_getEntryPtr_const')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! bool remove( std::string const& name, bool throwIfNotExists = true );
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterList_remove ( CT_Teuchos_ParameterList_ID_t selfID, char const name[], boolean throwIfNotExists );

  logical(c_bool) function Teuchos_ParameterList_remove ( selfID, name, throwIfNotExists ) &
        bind(C,name='Teuchos_ParameterList_remove')
    import :: c_bool ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    logical(c_bool)                    ,intent(in)   ,value              :: throwIfNotExists
  end function


  ! Original C++ prototype:
  ! ParameterList& sublist( const std::string& name, bool mustAlreadyExist = false ,const std::string& docString = "" );
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], boolean mustAlreadyExist, const char docString[] );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_sublist ( selfID, name, &
        mustAlreadyExist, docString ) bind(C,name='Teuchos_ParameterList_sublist')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char ,c_bool
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    logical(c_bool)                    ,intent(in)   ,value              :: mustAlreadyExist
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: docString
  end function


  ! Original C++ prototype:
  ! const ParameterList& sublist(const std::string& name) const;
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterList_sublist_existing ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterList_sublist_existing ( &
        selfID, name ) bind(C,name='Teuchos_ParameterList_sublist_existing')
    import :: FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! const std::string& name() const;
  ! CTrilinos prototype:
  ! const char * Teuchos_ParameterList_name_it ( CT_Teuchos_ParameterList_ID_t selfID );

  character(kind=c_char) function Teuchos_ParameterList_name_it ( selfID ) &
        bind(C,name='Teuchos_ParameterList_name_it')
    import :: kind=c_char ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! bool isParameter(const std::string& name) const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterList_isParameter ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  logical(c_bool) function Teuchos_ParameterList_isParameter ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_isParameter')
    import :: c_bool ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! bool isSublist(const std::string& name) const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterList_isSublist ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  logical(c_bool) function Teuchos_ParameterList_isSublist ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_isSublist')
    import :: c_bool ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> bool isType(const std::string& name) const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterList_isType_double ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  logical(c_bool) function Teuchos_ParameterList_isType_double ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_isType_double')
    import :: c_bool ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> bool isType(const std::string& name) const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterList_isType_int ( CT_Teuchos_ParameterList_ID_t selfID, const char name[] );

  logical(c_bool) function Teuchos_ParameterList_isType_int ( selfID, name ) &
        bind(C,name='Teuchos_ParameterList_isType_int')
    import :: c_bool ,FT_Teuchos_ParameterList_ID_t ,c_char
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
  end function


  ! Original C++ prototype:
  ! template<typename T> bool isType(const std::string& name, T* ptr) const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterList_isType_double_type ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], double * ptr );

  logical(c_bool) function Teuchos_ParameterList_isType_double_type ( selfID, name, ptr ) &
        bind(C,name='Teuchos_ParameterList_isType_double_type')
    import :: c_bool ,FT_Teuchos_ParameterList_ID_t ,c_char ,c_double
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    real(c_double)                                         ,dimension(*) :: ptr
  end function


  ! Original C++ prototype:
  ! template<typename T> bool isType(const std::string& name, T* ptr) const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterList_isType_int_type ( CT_Teuchos_ParameterList_ID_t selfID, const char name[], int * ptr );

  logical(c_bool) function Teuchos_ParameterList_isType_int_type ( selfID, name, ptr ) &
        bind(C,name='Teuchos_ParameterList_isType_int_type')
    import :: c_bool ,FT_Teuchos_ParameterList_ID_t ,c_char ,c_int
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)             ,intent(in)         ,dimension(*) :: name
    integer(c_int)                                         ,dimension(*) :: ptr
  end function


  ! Original C++ prototype:
  ! std::string currentParametersString() const;
  ! CTrilinos prototype:
  ! const char * Teuchos_ParameterList_currentParametersString ( CT_Teuchos_ParameterList_ID_t selfID );

  character(kind=c_char) function Teuchos_ParameterList_currentParametersString ( selfID ) &
        bind(C,name='Teuchos_ParameterList_currentParametersString')
    import :: kind=c_char ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! void validateParameters( ParameterList const& validParamList, int const depth = 1000, EValidateUsed const validateUsed = VALIDATE_USED_ENABLED, EValidateDefaults const validateDefaults = VALIDATE_DEFAULTS_ENABLED ) const;
  ! CTrilinos prototype:
  ! void Teuchos_ParameterList_validateParameters ( CT_Teuchos_ParameterList_ID_t selfID, CT_Teuchos_ParameterList_ID_t validParamListID, int const depth, const CT_EValidateUsed_E_t validateUsed, const CT_EValidateDefaults_E_t validateDefaults );

  subroutine Teuchos_ParameterList_validateParameters ( selfID, validParamListID, depth, &
        validateUsed, validateDefaults ) &
        bind(C,name='Teuchos_ParameterList_validateParameters')
    import :: FT_Teuchos_ParameterList_ID_t ,c_int ,FT_EValidateUsed_E_t , &
          FT_EValidateDefaults_E_t
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: validParamListID
    integer(c_int)                     ,intent(in)   ,value              :: depth
    integer(FT_EValidateUsed_E_t)      ,intent(in)   ,value              :: validateUsed
    integer(FT_EValidateDefaults_E_t)  ,intent(in)   ,value              :: validateDefaults
  end subroutine


  ! Original C++ prototype:
  ! void validateParametersAndSetDefaults( ParameterList const& validParamList, int const depth = 1000 );
  ! CTrilinos prototype:
  ! void Teuchos_ParameterList_validateParametersAndSetDefaults ( CT_Teuchos_ParameterList_ID_t selfID, CT_Teuchos_ParameterList_ID_t validParamListID, int const depth );

  subroutine Teuchos_ParameterList_validateParametersAndSetDefaults ( selfID, &
        validParamListID, depth ) &
        bind(C,name='Teuchos_ParameterList_validateParametersAndSetDefaults')
    import :: FT_Teuchos_ParameterList_ID_t ,c_int
    
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: validParamListID
    integer(c_int)                     ,intent(in)   ,value              :: depth
  end subroutine


  ! _________________ ParameterEntry interface bodies _________________


  ! Original C++ prototype:
  ! ParameterEntry();
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Create (  );

  type(FT_Teuchos_ParameterEntry_ID_t) function Teuchos_ParameterEntry_Create (  ) &
        bind(C,name='Teuchos_ParameterEntry_Create')
    import :: FT_Teuchos_ParameterEntry_ID_t
    
  end function


  ! Original C++ prototype:
  ! ParameterEntry(const ParameterEntry& source);
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterEntry_ID_t Teuchos_ParameterEntry_Duplicate ( CT_Teuchos_ParameterEntry_ID_t sourceID );

  type(FT_Teuchos_ParameterEntry_ID_t) function Teuchos_ParameterEntry_Duplicate ( sourceID ) &
        bind(C,name='Teuchos_ParameterEntry_Duplicate')
    import :: FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: sourceID
  end function


  ! Original C++ prototype:
  ! ~ParameterEntry();
  ! CTrilinos prototype:
  ! void Teuchos_ParameterEntry_Destroy ( CT_Teuchos_ParameterEntry_ID_t * selfID );

  subroutine Teuchos_ParameterEntry_Destroy ( selfID ) &
        bind(C,name='Teuchos_ParameterEntry_Destroy')
    import :: FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t)                                  :: selfID
  end subroutine


  ! Original C++ prototype:
  ! ParameterEntry& operator=(const ParameterEntry& source);
  ! CTrilinos prototype:
  ! void Teuchos_ParameterEntry_Assign ( CT_Teuchos_ParameterEntry_ID_t selfID, CT_Teuchos_ParameterEntry_ID_t sourceID );

  subroutine Teuchos_ParameterEntry_Assign ( selfID, sourceID ) &
        bind(C,name='Teuchos_ParameterEntry_Assign')
    import :: FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: sourceID
  end subroutine


  ! Original C++ prototype:
  ! void setAnyValue( const any &value, bool isDefault = false );
  ! CTrilinos prototype:
  ! void Teuchos_ParameterEntry_setAnyValue ( CT_Teuchos_ParameterEntry_ID_t selfID, CT_Teuchos_any_ID_t valueID, boolean isDefault );

  subroutine Teuchos_ParameterEntry_setAnyValue ( selfID, valueID, isDefault ) &
        bind(C,name='Teuchos_ParameterEntry_setAnyValue')
    import :: FT_Teuchos_ParameterEntry_ID_t ,FT_Teuchos_any_ID_t ,c_bool
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_any_ID_t)           ,intent(in)   ,value              :: valueID
    logical(c_bool)                     ,intent(in)   ,value              :: isDefault
  end subroutine


  ! Original C++ prototype:
  ! void setDocString(const std::string &docString);
  ! CTrilinos prototype:
  ! void Teuchos_ParameterEntry_setDocString ( CT_Teuchos_ParameterEntry_ID_t selfID, const char docString[] );

  subroutine Teuchos_ParameterEntry_setDocString ( selfID, docString ) &
        bind(C,name='Teuchos_ParameterEntry_setDocString')
    import :: FT_Teuchos_ParameterEntry_ID_t ,c_char
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    character(kind=c_char)              ,intent(in)         ,dimension(*) :: docString
  end subroutine


  ! Original C++ prototype:
  ! ParameterList& setList( bool isDefault = false, const std::string &docString = "" );
  ! CTrilinos prototype:
  ! CT_Teuchos_ParameterList_ID_t Teuchos_ParameterEntry_setList ( CT_Teuchos_ParameterEntry_ID_t selfID, boolean isDefault, const char docString[] );

  type(FT_Teuchos_ParameterList_ID_t) function Teuchos_ParameterEntry_setList ( selfID, &
        isDefault, docString ) bind(C,name='Teuchos_ParameterEntry_setList')
    import :: FT_Teuchos_ParameterList_ID_t ,FT_Teuchos_ParameterEntry_ID_t ,c_bool , &
          c_char
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    logical(c_bool)                     ,intent(in)   ,value              :: isDefault
    character(kind=c_char)              ,intent(in)         ,dimension(*) :: docString
  end function


  ! Original C++ prototype:
  ! template<typename T> inline T& getValue(T *ptr) const;
  ! CTrilinos prototype:
  ! double Teuchos_ParameterEntry_getValue_double ( CT_Teuchos_ParameterEntry_ID_t selfID, double * ptr );

  real(c_double) function Teuchos_ParameterEntry_getValue_double ( selfID, ptr ) &
        bind(C,name='Teuchos_ParameterEntry_getValue_double')
    import :: c_double ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    real(c_double)                                          ,dimension(*) :: ptr
  end function


  ! Original C++ prototype:
  ! template<typename T> inline T& getValue(T *ptr) const;
  ! CTrilinos prototype:
  ! int Teuchos_ParameterEntry_getValue_int ( CT_Teuchos_ParameterEntry_ID_t selfID, int * ptr );

  integer(c_int) function Teuchos_ParameterEntry_getValue_int ( selfID, ptr ) &
        bind(C,name='Teuchos_ParameterEntry_getValue_int')
    import :: c_int ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    integer(c_int)                                          ,dimension(*) :: ptr
  end function


  ! Original C++ prototype:
  ! inline any& getAny(bool activeQry = true);
  ! CTrilinos prototype:
  ! CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny ( CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );

  type(FT_Teuchos_any_ID_t) function Teuchos_ParameterEntry_getAny ( selfID, activeQry ) &
        bind(C,name='Teuchos_ParameterEntry_getAny')
    import :: FT_Teuchos_any_ID_t ,FT_Teuchos_ParameterEntry_ID_t ,c_bool
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    logical(c_bool)                     ,intent(in)   ,value              :: activeQry
  end function


  ! Original C++ prototype:
  ! inline const any& getAny(bool activeQry = true) const;
  ! CTrilinos prototype:
  ! CT_Teuchos_any_ID_t Teuchos_ParameterEntry_getAny_const ( CT_Teuchos_ParameterEntry_ID_t selfID, boolean activeQry );

  type(FT_Teuchos_any_ID_t) function Teuchos_ParameterEntry_getAny_const ( selfID, &
        activeQry ) bind(C,name='Teuchos_ParameterEntry_getAny_const')
    import :: FT_Teuchos_any_ID_t ,FT_Teuchos_ParameterEntry_ID_t ,c_bool
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
    logical(c_bool)                     ,intent(in)   ,value              :: activeQry
  end function


  ! Original C++ prototype:
  ! inline bool isUsed() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterEntry_isUsed ( CT_Teuchos_ParameterEntry_ID_t selfID );

  logical(c_bool) function Teuchos_ParameterEntry_isUsed ( selfID ) &
        bind(C,name='Teuchos_ParameterEntry_isUsed')
    import :: c_bool ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! bool isList() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterEntry_isList ( CT_Teuchos_ParameterEntry_ID_t selfID );

  logical(c_bool) function Teuchos_ParameterEntry_isList ( selfID ) &
        bind(C,name='Teuchos_ParameterEntry_isList')
    import :: c_bool ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! template <typename T> inline bool isType() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterEntry_isType_double ( CT_Teuchos_ParameterEntry_ID_t selfID );

  logical(c_bool) function Teuchos_ParameterEntry_isType_double ( selfID ) &
        bind(C,name='Teuchos_ParameterEntry_isType_double')
    import :: c_bool ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! template <typename T> inline bool isType() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterEntry_isType_int ( CT_Teuchos_ParameterEntry_ID_t selfID );

  logical(c_bool) function Teuchos_ParameterEntry_isType_int ( selfID ) &
        bind(C,name='Teuchos_ParameterEntry_isType_int')
    import :: c_bool ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! inline bool isDefault() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_ParameterEntry_isDefault ( CT_Teuchos_ParameterEntry_ID_t selfID );

  logical(c_bool) function Teuchos_ParameterEntry_isDefault ( selfID ) &
        bind(C,name='Teuchos_ParameterEntry_isDefault')
    import :: c_bool ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! inline std::string docString() const;
  ! CTrilinos prototype:
  ! const char * Teuchos_ParameterEntry_docString ( CT_Teuchos_ParameterEntry_ID_t selfID );

  character(kind=c_char) function Teuchos_ParameterEntry_docString ( selfID ) &
        bind(C,name='Teuchos_ParameterEntry_docString')
    import :: kind=c_char ,FT_Teuchos_ParameterEntry_ID_t
    
    type(FT_Teuchos_ParameterEntry_ID_t),intent(in)   ,value              :: selfID
  end function


  ! _________________ any interface bodies _________________


  ! Original C++ prototype:
  ! any();
  ! CTrilinos prototype:
  ! CT_Teuchos_any_ID_t Teuchos_any_Create (  );

  type(FT_Teuchos_any_ID_t) function Teuchos_any_Create (  ) &
        bind(C,name='Teuchos_any_Create')
    import :: FT_Teuchos_any_ID_t
    
  end function


  ! Original C++ prototype:
  ! any(const any & other);
  ! CTrilinos prototype:
  ! CT_Teuchos_any_ID_t Teuchos_any_Duplicate ( CT_Teuchos_any_ID_t otherID );

  type(FT_Teuchos_any_ID_t) function Teuchos_any_Duplicate ( otherID ) &
        bind(C,name='Teuchos_any_Duplicate')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: otherID
  end function


  ! Original C++ prototype:
  ! ~any();
  ! CTrilinos prototype:
  ! void Teuchos_any_Destroy ( CT_Teuchos_any_ID_t * selfID );

  subroutine Teuchos_any_Destroy ( selfID ) bind(C,name='Teuchos_any_Destroy')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)                                     :: selfID
  end subroutine


  ! Original C++ prototype:
  ! any & swap(any & rhs);
  ! CTrilinos prototype:
  ! CT_Teuchos_any_ID_t Teuchos_any_swap ( CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

  type(FT_Teuchos_any_ID_t) function Teuchos_any_swap ( selfID, rhsID ) &
        bind(C,name='Teuchos_any_swap')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: rhsID
  end function


  ! Original C++ prototype:
  ! any & operator=(const any & rhs);
  ! CTrilinos prototype:
  ! void Teuchos_any_Assign ( CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t rhsID );

  subroutine Teuchos_any_Assign ( selfID, rhsID ) bind(C,name='Teuchos_any_Assign')
    import :: FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: rhsID
  end subroutine


  ! Original C++ prototype:
  ! bool empty() const;
  ! CTrilinos prototype:
  ! boolean Teuchos_any_empty ( CT_Teuchos_any_ID_t selfID );

  logical(c_bool) function Teuchos_any_empty ( selfID ) bind(C,name='Teuchos_any_empty')
    import :: c_bool ,FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! std::string typeName() const;
  ! CTrilinos prototype:
  ! const char * Teuchos_any_typeName ( CT_Teuchos_any_ID_t selfID );

  character(kind=c_char) function Teuchos_any_typeName ( selfID ) &
        bind(C,name='Teuchos_any_typeName')
    import :: kind=c_char ,FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
  end function


  ! Original C++ prototype:
  ! bool same( const any &other ) const;
  ! CTrilinos prototype:
  ! boolean Teuchos_any_same ( CT_Teuchos_any_ID_t selfID, CT_Teuchos_any_ID_t otherID );

  logical(c_bool) function Teuchos_any_same ( selfID, otherID ) &
        bind(C,name='Teuchos_any_same')
    import :: c_bool ,FT_Teuchos_any_ID_t
    
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: selfID
    type(FT_Teuchos_any_ID_t)   ,intent(in)   ,value              :: otherID
  end function


  end interface
end module forteuchos

