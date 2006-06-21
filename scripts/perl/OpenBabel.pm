# This file was created automatically by SWIG 1.3.29.
# Don't modify this file, modify the SWIG interface instead.
package Chemistry::OpenBabel;
require Exporter;
require DynaLoader;
@ISA = qw(Exporter DynaLoader);
package Chemistry::OpenBabelc;
bootstrap Chemistry::OpenBabel;
package Chemistry::OpenBabel;
@EXPORT = qw( ); sub dl_load_flags { 0x01 }

# ---------- BASE METHODS -------------

package Chemistry::OpenBabel;

sub TIEHASH {
    my ($classname,$obj) = @_;
    return bless $obj, $classname;
}

sub CLEAR { }

sub FIRSTKEY { }

sub NEXTKEY { }

sub FETCH {
    my ($self,$field) = @_;
    my $member_func = "swig_${field}_get";
    $self->$member_func();
}

sub STORE {
    my ($self,$field,$newval) = @_;
    my $member_func = "swig_${field}_set";
    $self->$member_func($newval);
}

sub this {
    my $ptr = shift;
    return tied(%$ptr);
}


# ------- FUNCTION WRAPPERS --------

package Chemistry::OpenBabel;

*OpenDatafile = *Chemistry::OpenBabelc::OpenDatafile;
*DoubleMultiply = *Chemistry::OpenBabelc::DoubleMultiply;
*DoubleAdd = *Chemistry::OpenBabelc::DoubleAdd;
*DoubleModulus = *Chemistry::OpenBabelc::DoubleModulus;
*Point2Plane = *Chemistry::OpenBabelc::Point2Plane;
*Trim = *Chemistry::OpenBabelc::Trim;
*tokenize = *Chemistry::OpenBabelc::tokenize;
*ThrowError = *Chemistry::OpenBabelc::ThrowError;
*CartesianToInternal = *Chemistry::OpenBabelc::CartesianToInternal;
*InternalToCartesian = *Chemistry::OpenBabelc::InternalToCartesian;
*NewExtension = *Chemistry::OpenBabelc::NewExtension;
*get_rmat = *Chemistry::OpenBabelc::get_rmat;
*ob_make_rmat = *Chemistry::OpenBabelc::ob_make_rmat;
*qtrfit = *Chemistry::OpenBabelc::qtrfit;
*superimpose = *Chemistry::OpenBabelc::superimpose;
*CompareRingSize = *Chemistry::OpenBabelc::CompareRingSize;
*SmartsLexReplace = *Chemistry::OpenBabelc::SmartsLexReplace;

############# Class : Chemistry::OpenBabel::vectorInt ##############

package Chemistry::OpenBabel::vectorInt;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorInt(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorInt_size;
*empty = *Chemistry::OpenBabelc::vectorInt_empty;
*clear = *Chemistry::OpenBabelc::vectorInt_clear;
*push = *Chemistry::OpenBabelc::vectorInt_push;
*pop = *Chemistry::OpenBabelc::vectorInt_pop;
*get = *Chemistry::OpenBabelc::vectorInt_get;
*set = *Chemistry::OpenBabelc::vectorInt_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorInt($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vvInt ##############

package Chemistry::OpenBabel::vvInt;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vvInt(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vvInt_size;
*empty = *Chemistry::OpenBabelc::vvInt_empty;
*clear = *Chemistry::OpenBabelc::vvInt_clear;
*push = *Chemistry::OpenBabelc::vvInt_push;
*pop = *Chemistry::OpenBabelc::vvInt_pop;
*get = *Chemistry::OpenBabelc::vvInt_get;
*set = *Chemistry::OpenBabelc::vvInt_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vvInt($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorDouble ##############

package Chemistry::OpenBabel::vectorDouble;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorDouble(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorDouble_size;
*empty = *Chemistry::OpenBabelc::vectorDouble_empty;
*clear = *Chemistry::OpenBabelc::vectorDouble_clear;
*push = *Chemistry::OpenBabelc::vectorDouble_push;
*pop = *Chemistry::OpenBabelc::vectorDouble_pop;
*get = *Chemistry::OpenBabelc::vectorDouble_get;
*set = *Chemistry::OpenBabelc::vectorDouble_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorDouble($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vVector3 ##############

package Chemistry::OpenBabel::vVector3;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vVector3(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vVector3_size;
*empty = *Chemistry::OpenBabelc::vVector3_empty;
*clear = *Chemistry::OpenBabelc::vVector3_clear;
*push = *Chemistry::OpenBabelc::vVector3_push;
*pop = *Chemistry::OpenBabelc::vVector3_pop;
*get = *Chemistry::OpenBabelc::vVector3_get;
*set = *Chemistry::OpenBabelc::vVector3_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vVector3($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorMol ##############

package Chemistry::OpenBabel::vectorMol;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorMol(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorMol_size;
*empty = *Chemistry::OpenBabelc::vectorMol_empty;
*clear = *Chemistry::OpenBabelc::vectorMol_clear;
*push = *Chemistry::OpenBabelc::vectorMol_push;
*pop = *Chemistry::OpenBabelc::vectorMol_pop;
*get = *Chemistry::OpenBabelc::vectorMol_get;
*set = *Chemistry::OpenBabelc::vectorMol_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorMol($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorBond ##############

package Chemistry::OpenBabel::vectorBond;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorBond(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorBond_size;
*empty = *Chemistry::OpenBabelc::vectorBond_empty;
*clear = *Chemistry::OpenBabelc::vectorBond_clear;
*push = *Chemistry::OpenBabelc::vectorBond_push;
*pop = *Chemistry::OpenBabelc::vectorBond_pop;
*get = *Chemistry::OpenBabelc::vectorBond_get;
*set = *Chemistry::OpenBabelc::vectorBond_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorBond($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorResidue ##############

package Chemistry::OpenBabel::vectorResidue;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorResidue(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorResidue_size;
*empty = *Chemistry::OpenBabelc::vectorResidue_empty;
*clear = *Chemistry::OpenBabelc::vectorResidue_clear;
*push = *Chemistry::OpenBabelc::vectorResidue_push;
*pop = *Chemistry::OpenBabelc::vectorResidue_pop;
*get = *Chemistry::OpenBabelc::vectorResidue_get;
*set = *Chemistry::OpenBabelc::vectorResidue_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorResidue($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorRing ##############

package Chemistry::OpenBabel::vectorRing;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorRing(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorRing_size;
*empty = *Chemistry::OpenBabelc::vectorRing_empty;
*clear = *Chemistry::OpenBabelc::vectorRing_clear;
*push = *Chemistry::OpenBabelc::vectorRing_push;
*pop = *Chemistry::OpenBabelc::vectorRing_pop;
*get = *Chemistry::OpenBabelc::vectorRing_get;
*set = *Chemistry::OpenBabelc::vectorRing_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorRing($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vectorData ##############

package Chemistry::OpenBabel::vectorData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vectorData(@_);
    bless $self, $pkg if defined($self);
}

*size = *Chemistry::OpenBabelc::vectorData_size;
*empty = *Chemistry::OpenBabelc::vectorData_empty;
*clear = *Chemistry::OpenBabelc::vectorData_clear;
*push = *Chemistry::OpenBabelc::vectorData_push;
*pop = *Chemistry::OpenBabelc::vectorData_pop;
*get = *Chemistry::OpenBabelc::vectorData_get;
*set = *Chemistry::OpenBabelc::vectorData_set;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vectorData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBGlobalDataBase ##############

package Chemistry::OpenBabel::OBGlobalDataBase;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBGlobalDataBase(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBGlobalDataBase($self);
        delete $OWNER{$self};
    }
}

*Init = *Chemistry::OpenBabelc::OBGlobalDataBase_Init;
*GetSize = *Chemistry::OpenBabelc::OBGlobalDataBase_GetSize;
*SetReadDirectory = *Chemistry::OpenBabelc::OBGlobalDataBase_SetReadDirectory;
*SetEnvironmentVariable = *Chemistry::OpenBabelc::OBGlobalDataBase_SetEnvironmentVariable;
*ParseLine = *Chemistry::OpenBabelc::OBGlobalDataBase_ParseLine;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBElement ##############

package Chemistry::OpenBabel::OBElement;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBElement(@_);
    bless $self, $pkg if defined($self);
}

*GetAtomicNum = *Chemistry::OpenBabelc::OBElement_GetAtomicNum;
*GetSymbol = *Chemistry::OpenBabelc::OBElement_GetSymbol;
*GetCovalentRad = *Chemistry::OpenBabelc::OBElement_GetCovalentRad;
*GetVdwRad = *Chemistry::OpenBabelc::OBElement_GetVdwRad;
*GetMass = *Chemistry::OpenBabelc::OBElement_GetMass;
*GetMaxBonds = *Chemistry::OpenBabelc::OBElement_GetMaxBonds;
*GetElectroNeg = *Chemistry::OpenBabelc::OBElement_GetElectroNeg;
*GetIonization = *Chemistry::OpenBabelc::OBElement_GetIonization;
*GetElectronAffinity = *Chemistry::OpenBabelc::OBElement_GetElectronAffinity;
*GetName = *Chemistry::OpenBabelc::OBElement_GetName;
*GetRed = *Chemistry::OpenBabelc::OBElement_GetRed;
*GetGreen = *Chemistry::OpenBabelc::OBElement_GetGreen;
*GetBlue = *Chemistry::OpenBabelc::OBElement_GetBlue;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBElement($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBElementTable ##############

package Chemistry::OpenBabel::OBElementTable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGlobalDataBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBElementTable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBElementTable($self);
        delete $OWNER{$self};
    }
}

*ParseLine = *Chemistry::OpenBabelc::OBElementTable_ParseLine;
*GetNumberOfElements = *Chemistry::OpenBabelc::OBElementTable_GetNumberOfElements;
*GetSize = *Chemistry::OpenBabelc::OBElementTable_GetSize;
*GetAtomicNum = *Chemistry::OpenBabelc::OBElementTable_GetAtomicNum;
*GetSymbol = *Chemistry::OpenBabelc::OBElementTable_GetSymbol;
*GetVdwRad = *Chemistry::OpenBabelc::OBElementTable_GetVdwRad;
*GetCovalentRad = *Chemistry::OpenBabelc::OBElementTable_GetCovalentRad;
*GetMass = *Chemistry::OpenBabelc::OBElementTable_GetMass;
*CorrectedBondRad = *Chemistry::OpenBabelc::OBElementTable_CorrectedBondRad;
*CorrectedVdwRad = *Chemistry::OpenBabelc::OBElementTable_CorrectedVdwRad;
*GetMaxBonds = *Chemistry::OpenBabelc::OBElementTable_GetMaxBonds;
*GetElectroNeg = *Chemistry::OpenBabelc::OBElementTable_GetElectroNeg;
*GetIonization = *Chemistry::OpenBabelc::OBElementTable_GetIonization;
*GetElectronAffinity = *Chemistry::OpenBabelc::OBElementTable_GetElectronAffinity;
*GetRGB = *Chemistry::OpenBabelc::OBElementTable_GetRGB;
*GetName = *Chemistry::OpenBabelc::OBElementTable_GetName;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBIsotopeTable ##############

package Chemistry::OpenBabel::OBIsotopeTable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGlobalDataBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBIsotopeTable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBIsotopeTable($self);
        delete $OWNER{$self};
    }
}

*GetSize = *Chemistry::OpenBabelc::OBIsotopeTable_GetSize;
*ParseLine = *Chemistry::OpenBabelc::OBIsotopeTable_ParseLine;
*GetExactMass = *Chemistry::OpenBabelc::OBIsotopeTable_GetExactMass;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBTypeTable ##############

package Chemistry::OpenBabel::OBTypeTable;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGlobalDataBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBTypeTable(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBTypeTable($self);
        delete $OWNER{$self};
    }
}

*ParseLine = *Chemistry::OpenBabelc::OBTypeTable_ParseLine;
*GetSize = *Chemistry::OpenBabelc::OBTypeTable_GetSize;
*SetFromType = *Chemistry::OpenBabelc::OBTypeTable_SetFromType;
*SetToType = *Chemistry::OpenBabelc::OBTypeTable_SetToType;
*Translate = *Chemistry::OpenBabelc::OBTypeTable_Translate;
*GetFromType = *Chemistry::OpenBabelc::OBTypeTable_GetFromType;
*GetToType = *Chemistry::OpenBabelc::OBTypeTable_GetToType;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBResidueData ##############

package Chemistry::OpenBabel::OBResidueData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGlobalDataBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBResidueData(@_);
    bless $self, $pkg if defined($self);
}

*ParseLine = *Chemistry::OpenBabelc::OBResidueData_ParseLine;
*GetSize = *Chemistry::OpenBabelc::OBResidueData_GetSize;
*SetResName = *Chemistry::OpenBabelc::OBResidueData_SetResName;
*LookupBO = *Chemistry::OpenBabelc::OBResidueData_LookupBO;
*LookupType = *Chemistry::OpenBabelc::OBResidueData_LookupType;
*AssignBonds = *Chemistry::OpenBabelc::OBResidueData_AssignBonds;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBResidueData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBStopwatch ##############

package Chemistry::OpenBabel::OBStopwatch;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Start = *Chemistry::OpenBabelc::OBStopwatch_Start;
*Lap = *Chemistry::OpenBabelc::OBStopwatch_Lap;
*Elapsed = *Chemistry::OpenBabelc::OBStopwatch_Elapsed;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBStopwatch(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBStopwatch($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSqrtTbl ##############

package Chemistry::OpenBabel::OBSqrtTbl;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSqrtTbl(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSqrtTbl($self);
        delete $OWNER{$self};
    }
}

*Sqrt = *Chemistry::OpenBabelc::OBSqrtTbl_Sqrt;
*Init = *Chemistry::OpenBabelc::OBSqrtTbl_Init;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::DoubleType ##############

package Chemistry::OpenBabel::DoubleType;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig_hi_get = *Chemistry::OpenBabelc::DoubleType_hi_get;
*swig_hi_set = *Chemistry::OpenBabelc::DoubleType_hi_set;
*swig_lo_get = *Chemistry::OpenBabelc::DoubleType_lo_get;
*swig_lo_set = *Chemistry::OpenBabelc::DoubleType_lo_set;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_DoubleType(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_DoubleType($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRandom ##############

package Chemistry::OpenBabel::OBRandom;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRandom(@_);
    bless $self, $pkg if defined($self);
}

*Seed = *Chemistry::OpenBabelc::OBRandom_Seed;
*TimeSeed = *Chemistry::OpenBabelc::OBRandom_TimeSeed;
*NextInt = *Chemistry::OpenBabelc::OBRandom_NextInt;
*NextFloat = *Chemistry::OpenBabelc::OBRandom_NextFloat;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRandom($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::vector3 ##############

package Chemistry::OpenBabel::vector3;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_vector3(@_);
    bless $self, $pkg if defined($self);
}

*Set = *Chemistry::OpenBabelc::vector3_Set;
*SetX = *Chemistry::OpenBabelc::vector3_SetX;
*SetY = *Chemistry::OpenBabelc::vector3_SetY;
*SetZ = *Chemistry::OpenBabelc::vector3_SetZ;
*Get = *Chemistry::OpenBabelc::vector3_Get;
*randomUnitVector = *Chemistry::OpenBabelc::vector3_randomUnitVector;
*normalize = *Chemistry::OpenBabelc::vector3_normalize;
*length = *Chemistry::OpenBabelc::vector3_length;
*length_2 = *Chemistry::OpenBabelc::vector3_length_2;
*x = *Chemistry::OpenBabelc::vector3_x;
*y = *Chemistry::OpenBabelc::vector3_y;
*z = *Chemistry::OpenBabelc::vector3_z;
*distSq = *Chemistry::OpenBabelc::vector3_distSq;
*createOrthoVector = *Chemistry::OpenBabelc::vector3_createOrthoVector;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_vector3($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBGenericData ##############

package Chemistry::OpenBabel::OBGenericData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBGenericData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBGenericData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBGenericData($self);
        delete $OWNER{$self};
    }
}

*SetAttribute = *Chemistry::OpenBabelc::OBGenericData_SetAttribute;
*GetAttribute = *Chemistry::OpenBabelc::OBGenericData_GetAttribute;
*GetDataType = *Chemistry::OpenBabelc::OBGenericData_GetDataType;
*GetValue = *Chemistry::OpenBabelc::OBGenericData_GetValue;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBCommentData ##############

package Chemistry::OpenBabel::OBCommentData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBCommentData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBCommentData_Clone;
*SetData = *Chemistry::OpenBabelc::OBCommentData_SetData;
*GetData = *Chemistry::OpenBabelc::OBCommentData_GetData;
*GetValue = *Chemistry::OpenBabelc::OBCommentData_GetValue;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBCommentData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBExternalBond ##############

package Chemistry::OpenBabel::OBExternalBond;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBExternalBond(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBExternalBond($self);
        delete $OWNER{$self};
    }
}

*GetIdx = *Chemistry::OpenBabelc::OBExternalBond_GetIdx;
*GetAtom = *Chemistry::OpenBabelc::OBExternalBond_GetAtom;
*GetBond = *Chemistry::OpenBabelc::OBExternalBond_GetBond;
*SetIdx = *Chemistry::OpenBabelc::OBExternalBond_SetIdx;
*SetAtom = *Chemistry::OpenBabelc::OBExternalBond_SetAtom;
*SetBond = *Chemistry::OpenBabelc::OBExternalBond_SetBond;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBExternalBondData ##############

package Chemistry::OpenBabel::OBExternalBondData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBExternalBondData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBExternalBondData_Clone;
*SetData = *Chemistry::OpenBabelc::OBExternalBondData_SetData;
*GetData = *Chemistry::OpenBabelc::OBExternalBondData_GetData;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBExternalBondData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBPairData ##############

package Chemistry::OpenBabel::OBPairData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBPairData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBPairData_Clone;
*SetValue = *Chemistry::OpenBabelc::OBPairData_SetValue;
*GetValue = *Chemistry::OpenBabelc::OBPairData_GetValue;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBPairData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBVirtualBond ##############

package Chemistry::OpenBabel::OBVirtualBond;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Clone = *Chemistry::OpenBabelc::OBVirtualBond_Clone;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBVirtualBond(@_);
    bless $self, $pkg if defined($self);
}

*GetBgn = *Chemistry::OpenBabelc::OBVirtualBond_GetBgn;
*GetEnd = *Chemistry::OpenBabelc::OBVirtualBond_GetEnd;
*GetOrder = *Chemistry::OpenBabelc::OBVirtualBond_GetOrder;
*GetStereo = *Chemistry::OpenBabelc::OBVirtualBond_GetStereo;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBVirtualBond($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRingData ##############

package Chemistry::OpenBabel::OBRingData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRingData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBRingData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRingData($self);
        delete $OWNER{$self};
    }
}

*SetData = *Chemistry::OpenBabelc::OBRingData_SetData;
*PushBack = *Chemistry::OpenBabelc::OBRingData_PushBack;
*GetData = *Chemistry::OpenBabelc::OBRingData_GetData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBUnitCell ##############

package Chemistry::OpenBabel::OBUnitCell;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBUnitCell(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBUnitCell_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBUnitCell($self);
        delete $OWNER{$self};
    }
}

*SetData = *Chemistry::OpenBabelc::OBUnitCell_SetData;
*SetOffset = *Chemistry::OpenBabelc::OBUnitCell_SetOffset;
*SetSpaceGroup = *Chemistry::OpenBabelc::OBUnitCell_SetSpaceGroup;
*GetA = *Chemistry::OpenBabelc::OBUnitCell_GetA;
*GetB = *Chemistry::OpenBabelc::OBUnitCell_GetB;
*GetC = *Chemistry::OpenBabelc::OBUnitCell_GetC;
*GetAlpha = *Chemistry::OpenBabelc::OBUnitCell_GetAlpha;
*GetBeta = *Chemistry::OpenBabelc::OBUnitCell_GetBeta;
*GetGamma = *Chemistry::OpenBabelc::OBUnitCell_GetGamma;
*GetOffset = *Chemistry::OpenBabelc::OBUnitCell_GetOffset;
*GetSpaceGroup = *Chemistry::OpenBabelc::OBUnitCell_GetSpaceGroup;
*GetCellVectors = *Chemistry::OpenBabelc::OBUnitCell_GetCellVectors;
*GetCellMatrix = *Chemistry::OpenBabelc::OBUnitCell_GetCellMatrix;
*GetOrthoMatrix = *Chemistry::OpenBabelc::OBUnitCell_GetOrthoMatrix;
*GetFractionalMatrix = *Chemistry::OpenBabelc::OBUnitCell_GetFractionalMatrix;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBConformerData ##############

package Chemistry::OpenBabel::OBConformerData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBConformerData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBConformerData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBConformerData($self);
        delete $OWNER{$self};
    }
}

*SetDimension = *Chemistry::OpenBabelc::OBConformerData_SetDimension;
*SetEnergies = *Chemistry::OpenBabelc::OBConformerData_SetEnergies;
*SetForces = *Chemistry::OpenBabelc::OBConformerData_SetForces;
*SetVelocities = *Chemistry::OpenBabelc::OBConformerData_SetVelocities;
*SetDisplacements = *Chemistry::OpenBabelc::OBConformerData_SetDisplacements;
*SetData = *Chemistry::OpenBabelc::OBConformerData_SetData;
*GetDimension = *Chemistry::OpenBabelc::OBConformerData_GetDimension;
*GetEnergies = *Chemistry::OpenBabelc::OBConformerData_GetEnergies;
*GetForces = *Chemistry::OpenBabelc::OBConformerData_GetForces;
*GetVelocities = *Chemistry::OpenBabelc::OBConformerData_GetVelocities;
*GetDisplacements = *Chemistry::OpenBabelc::OBConformerData_GetDisplacements;
*GetData = *Chemistry::OpenBabelc::OBConformerData_GetData;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSymmetryData ##############

package Chemistry::OpenBabel::OBSymmetryData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSymmetryData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBSymmetryData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSymmetryData($self);
        delete $OWNER{$self};
    }
}

*SetData = *Chemistry::OpenBabelc::OBSymmetryData_SetData;
*SetPointGroup = *Chemistry::OpenBabelc::OBSymmetryData_SetPointGroup;
*SetSpaceGroup = *Chemistry::OpenBabelc::OBSymmetryData_SetSpaceGroup;
*GetPointGroup = *Chemistry::OpenBabelc::OBSymmetryData_GetPointGroup;
*GetSpaceGroup = *Chemistry::OpenBabelc::OBSymmetryData_GetSpaceGroup;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBTorsion ##############

package Chemistry::OpenBabel::OBTorsion;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBTorsion(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBTorsion($self);
        delete $OWNER{$self};
    }
}

*Clear = *Chemistry::OpenBabelc::OBTorsion_Clear;
*Empty = *Chemistry::OpenBabelc::OBTorsion_Empty;
*AddTorsion = *Chemistry::OpenBabelc::OBTorsion_AddTorsion;
*SetAngle = *Chemistry::OpenBabelc::OBTorsion_SetAngle;
*SetData = *Chemistry::OpenBabelc::OBTorsion_SetData;
*GetAngle = *Chemistry::OpenBabelc::OBTorsion_GetAngle;
*GetBondIdx = *Chemistry::OpenBabelc::OBTorsion_GetBondIdx;
*GetSize = *Chemistry::OpenBabelc::OBTorsion_GetSize;
*GetBC = *Chemistry::OpenBabelc::OBTorsion_GetBC;
*GetADs = *Chemistry::OpenBabelc::OBTorsion_GetADs;
*IsProtonRotor = *Chemistry::OpenBabelc::OBTorsion_IsProtonRotor;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBTorsionData ##############

package Chemistry::OpenBabel::OBTorsionData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Clone = *Chemistry::OpenBabelc::OBTorsionData_Clone;
*Clear = *Chemistry::OpenBabelc::OBTorsionData_Clear;
*GetData = *Chemistry::OpenBabelc::OBTorsionData_GetData;
*GetSize = *Chemistry::OpenBabelc::OBTorsionData_GetSize;
*SetData = *Chemistry::OpenBabelc::OBTorsionData_SetData;
*FillTorsionArray = *Chemistry::OpenBabelc::OBTorsionData_FillTorsionArray;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBTorsionData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBAngle ##############

package Chemistry::OpenBabel::OBAngle;
use overload
    "==" => sub { $_[0]->__eq__($_[1])},
    "fallback" => 1;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBAngle(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBAngle($self);
        delete $OWNER{$self};
    }
}

*__eq__ = *Chemistry::OpenBabelc::OBAngle___eq__;
*Clear = *Chemistry::OpenBabelc::OBAngle_Clear;
*GetAngle = *Chemistry::OpenBabelc::OBAngle_GetAngle;
*SetAngle = *Chemistry::OpenBabelc::OBAngle_SetAngle;
*SetAtoms = *Chemistry::OpenBabelc::OBAngle_SetAtoms;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBAngleData ##############

package Chemistry::OpenBabel::OBAngleData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*Clone = *Chemistry::OpenBabelc::OBAngleData_Clone;
*Clear = *Chemistry::OpenBabelc::OBAngleData_Clear;
*FillAngleArray = *Chemistry::OpenBabelc::OBAngleData_FillAngleArray;
*SetData = *Chemistry::OpenBabelc::OBAngleData_SetData;
*GetSize = *Chemistry::OpenBabelc::OBAngleData_GetSize;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBAngleData($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBChiralData ##############

package Chemistry::OpenBabel::OBChiralData;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*GetAtom4Refs = *Chemistry::OpenBabelc::OBChiralData_GetAtom4Refs;
*GetAtomRef = *Chemistry::OpenBabelc::OBChiralData_GetAtomRef;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBChiralData(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBChiralData_Clone;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBChiralData($self);
        delete $OWNER{$self};
    }
}

*Clear = *Chemistry::OpenBabelc::OBChiralData_Clear;
*SetAtom4Refs = *Chemistry::OpenBabelc::OBChiralData_SetAtom4Refs;
*AddAtomRef = *Chemistry::OpenBabelc::OBChiralData_AddAtomRef;
*GetSize = *Chemistry::OpenBabelc::OBChiralData_GetSize;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSerialNums ##############

package Chemistry::OpenBabel::OBSerialNums;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGenericData Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSerialNums(@_);
    bless $self, $pkg if defined($self);
}

*Clone = *Chemistry::OpenBabelc::OBSerialNums_Clone;
*GetData = *Chemistry::OpenBabelc::OBSerialNums_GetData;
*SetData = *Chemistry::OpenBabelc::OBSerialNums_SetData;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSerialNums($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBBase ##############

package Chemistry::OpenBabel::OBBase;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBBase($self);
        delete $OWNER{$self};
    }
}

*DoTransformations = *Chemistry::OpenBabelc::OBBase_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBBase_ClassDescription;
*HasData = *Chemistry::OpenBabelc::OBBase_HasData;
*DeleteData = *Chemistry::OpenBabelc::OBBase_DeleteData;
*SetData = *Chemistry::OpenBabelc::OBBase_SetData;
*DataSize = *Chemistry::OpenBabelc::OBBase_DataSize;
*GetData = *Chemistry::OpenBabelc::OBBase_GetData;
*BeginData = *Chemistry::OpenBabelc::OBBase_BeginData;
*EndData = *Chemistry::OpenBabelc::OBBase_EndData;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBBase(@_);
    bless $self, $pkg if defined($self);
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBNodeBase ##############

package Chemistry::OpenBabel::OBNodeBase;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig_Visit_get = *Chemistry::OpenBabelc::OBNodeBase_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBNodeBase_Visit_set;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBNodeBase(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBNodeBase($self);
        delete $OWNER{$self};
    }
}

*GetIdx = *Chemistry::OpenBabelc::OBNodeBase_GetIdx;
*SetIdx = *Chemistry::OpenBabelc::OBNodeBase_SetIdx;
*GetParent = *Chemistry::OpenBabelc::OBNodeBase_GetParent;
*SetParent = *Chemistry::OpenBabelc::OBNodeBase_SetParent;
*AddEdge = *Chemistry::OpenBabelc::OBNodeBase_AddEdge;
*GetValence = *Chemistry::OpenBabelc::OBNodeBase_GetValence;
*IsConnected = *Chemistry::OpenBabelc::OBNodeBase_IsConnected;
*Error = *Chemistry::OpenBabelc::OBNodeBase_Error;
*GetFormalCharge = *Chemistry::OpenBabelc::OBNodeBase_GetFormalCharge;
*ExplicitHydrogenCount = *Chemistry::OpenBabelc::OBNodeBase_ExplicitHydrogenCount;
*ImplicitHydrogenCount = *Chemistry::OpenBabelc::OBNodeBase_ImplicitHydrogenCount;
*GetImplicitValence = *Chemistry::OpenBabelc::OBNodeBase_GetImplicitValence;
*GetHvyValence = *Chemistry::OpenBabelc::OBNodeBase_GetHvyValence;
*KBOSum = *Chemistry::OpenBabelc::OBNodeBase_KBOSum;
*GetHyb = *Chemistry::OpenBabelc::OBNodeBase_GetHyb;
*MemberOfRingCount = *Chemistry::OpenBabelc::OBNodeBase_MemberOfRingCount;
*GetAtomicNum = *Chemistry::OpenBabelc::OBNodeBase_GetAtomicNum;
*SetMatch = *Chemistry::OpenBabelc::OBNodeBase_SetMatch;
*SetAromatic = *Chemistry::OpenBabelc::OBNodeBase_SetAromatic;
*IsInRingSize = *Chemistry::OpenBabelc::OBNodeBase_IsInRingSize;
*IsAromatic = *Chemistry::OpenBabelc::OBNodeBase_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBNodeBase_IsInRing;
*Eval = *Chemistry::OpenBabelc::OBNodeBase_Eval;
*GetMatch = *Chemistry::OpenBabelc::OBNodeBase_GetMatch;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBEdgeBase ##############

package Chemistry::OpenBabel::OBEdgeBase;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig_Visit_get = *Chemistry::OpenBabelc::OBEdgeBase_Visit_get;
*swig_Visit_set = *Chemistry::OpenBabelc::OBEdgeBase_Visit_set;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBEdgeBase(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBEdgeBase($self);
        delete $OWNER{$self};
    }
}

*GetParent = *Chemistry::OpenBabelc::OBEdgeBase_GetParent;
*SetParent = *Chemistry::OpenBabelc::OBEdgeBase_SetParent;
*GetIdx = *Chemistry::OpenBabelc::OBEdgeBase_GetIdx;
*SetIdx = *Chemistry::OpenBabelc::OBEdgeBase_SetIdx;
*SetBgn = *Chemistry::OpenBabelc::OBEdgeBase_SetBgn;
*SetEnd = *Chemistry::OpenBabelc::OBEdgeBase_SetEnd;
*SwapEnds = *Chemistry::OpenBabelc::OBEdgeBase_SwapEnds;
*GetBgn = *Chemistry::OpenBabelc::OBEdgeBase_GetBgn;
*GetEnd = *Chemistry::OpenBabelc::OBEdgeBase_GetEnd;
*Error = *Chemistry::OpenBabelc::OBEdgeBase_Error;
*SetClosure = *Chemistry::OpenBabelc::OBEdgeBase_SetClosure;
*IsAromatic = *Chemistry::OpenBabelc::OBEdgeBase_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBEdgeBase_IsInRing;
*IsClosure = *Chemistry::OpenBabelc::OBEdgeBase_IsClosure;
*Eval = *Chemistry::OpenBabelc::OBEdgeBase_Eval;
*GetBO = *Chemistry::OpenBabelc::OBEdgeBase_GetBO;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBGraphBase ##############

package Chemistry::OpenBabel::OBGraphBase;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBGraphBase(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBGraphBase($self);
        delete $OWNER{$self};
    }
}

*NumNodes = *Chemistry::OpenBabelc::OBGraphBase_NumNodes;
*NumEdges = *Chemistry::OpenBabelc::OBGraphBase_NumEdges;
*ResetVisitFlags = *Chemistry::OpenBabelc::OBGraphBase_ResetVisitFlags;
*SetVisitLock = *Chemistry::OpenBabelc::OBGraphBase_SetVisitLock;
*GetVisitLock = *Chemistry::OpenBabelc::OBGraphBase_GetVisitLock;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBFormat ##############

package Chemistry::OpenBabel::OBFormat;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*ReadMolecule = *Chemistry::OpenBabelc::OBFormat_ReadMolecule;
*ReadChemObject = *Chemistry::OpenBabelc::OBFormat_ReadChemObject;
*WriteMolecule = *Chemistry::OpenBabelc::OBFormat_WriteMolecule;
*WriteChemObject = *Chemistry::OpenBabelc::OBFormat_WriteChemObject;
*Description = *Chemistry::OpenBabelc::OBFormat_Description;
*TargetClassDescription = *Chemistry::OpenBabelc::OBFormat_TargetClassDescription;
*GetType = *Chemistry::OpenBabelc::OBFormat_GetType;
*SpecificationURL = *Chemistry::OpenBabelc::OBFormat_SpecificationURL;
*GetMIMEType = *Chemistry::OpenBabelc::OBFormat_GetMIMEType;
*Flags = *Chemistry::OpenBabelc::OBFormat_Flags;
*SkipObjects = *Chemistry::OpenBabelc::OBFormat_SkipObjects;
*MakeNewInstance = *Chemistry::OpenBabelc::OBFormat_MakeNewInstance;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBFormat($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::CharPtrLess ##############

package Chemistry::OpenBabel::CharPtrLess;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*__call__ = *Chemistry::OpenBabelc::CharPtrLess___call__;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_CharPtrLess(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_CharPtrLess($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBConversion ##############

package Chemistry::OpenBabel::OBConversion;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBConversion(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBConversion($self);
        delete $OWNER{$self};
    }
}

*RegisterFormat = *Chemistry::OpenBabelc::OBConversion_RegisterFormat;
*FindFormat = *Chemistry::OpenBabelc::OBConversion_FindFormat;
*FormatFromExt = *Chemistry::OpenBabelc::OBConversion_FormatFromExt;
*FormatFromMIME = *Chemistry::OpenBabelc::OBConversion_FormatFromMIME;
*GetNextFormat = *Chemistry::OpenBabelc::OBConversion_GetNextFormat;
*Description = *Chemistry::OpenBabelc::OBConversion_Description;
*GetInStream = *Chemistry::OpenBabelc::OBConversion_GetInStream;
*GetOutStream = *Chemistry::OpenBabelc::OBConversion_GetOutStream;
*SetInStream = *Chemistry::OpenBabelc::OBConversion_SetInStream;
*SetOutStream = *Chemistry::OpenBabelc::OBConversion_SetOutStream;
*SetInAndOutFormats = *Chemistry::OpenBabelc::OBConversion_SetInAndOutFormats;
*SetInFormat = *Chemistry::OpenBabelc::OBConversion_SetInFormat;
*SetOutFormat = *Chemistry::OpenBabelc::OBConversion_SetOutFormat;
*GetInFormat = *Chemistry::OpenBabelc::OBConversion_GetInFormat;
*GetOutFormat = *Chemistry::OpenBabelc::OBConversion_GetOutFormat;
*GetInFilename = *Chemistry::OpenBabelc::OBConversion_GetInFilename;
*GetInPos = *Chemistry::OpenBabelc::OBConversion_GetInPos;
*GetInLen = *Chemistry::OpenBabelc::OBConversion_GetInLen;
*GetTitle = *Chemistry::OpenBabelc::OBConversion_GetTitle;
*GetAuxConv = *Chemistry::OpenBabelc::OBConversion_GetAuxConv;
*SetAuxConv = *Chemistry::OpenBabelc::OBConversion_SetAuxConv;
*INOPTIONS = *Chemistry::OpenBabelc::OBConversion_INOPTIONS;
*OUTOPTIONS = *Chemistry::OpenBabelc::OBConversion_OUTOPTIONS;
*GENOPTIONS = *Chemistry::OpenBabelc::OBConversion_GENOPTIONS;
*IsOption = *Chemistry::OpenBabelc::OBConversion_IsOption;
*GetOptions = *Chemistry::OpenBabelc::OBConversion_GetOptions;
*AddOption = *Chemistry::OpenBabelc::OBConversion_AddOption;
*RemoveOption = *Chemistry::OpenBabelc::OBConversion_RemoveOption;
*SetOptions = *Chemistry::OpenBabelc::OBConversion_SetOptions;
*RegisterOptionParam = *Chemistry::OpenBabelc::OBConversion_RegisterOptionParam;
*GetOptionParams = *Chemistry::OpenBabelc::OBConversion_GetOptionParams;
*Convert = *Chemistry::OpenBabelc::OBConversion_Convert;
*FullConvert = *Chemistry::OpenBabelc::OBConversion_FullConvert;
*AddChemObject = *Chemistry::OpenBabelc::OBConversion_AddChemObject;
*GetChemObject = *Chemistry::OpenBabelc::OBConversion_GetChemObject;
*IsLast = *Chemistry::OpenBabelc::OBConversion_IsLast;
*IsFirstInput = *Chemistry::OpenBabelc::OBConversion_IsFirstInput;
*GetOutputIndex = *Chemistry::OpenBabelc::OBConversion_GetOutputIndex;
*SetOutputIndex = *Chemistry::OpenBabelc::OBConversion_SetOutputIndex;
*SetMoreFilesToCome = *Chemistry::OpenBabelc::OBConversion_SetMoreFilesToCome;
*SetOneObjectOnly = *Chemistry::OpenBabelc::OBConversion_SetOneObjectOnly;
*GetDefaultFormat = *Chemistry::OpenBabelc::OBConversion_GetDefaultFormat;
*Write = *Chemistry::OpenBabelc::OBConversion_Write;
*WriteString = *Chemistry::OpenBabelc::OBConversion_WriteString;
*WriteFile = *Chemistry::OpenBabelc::OBConversion_WriteFile;
*Read = *Chemistry::OpenBabelc::OBConversion_Read;
*ReadString = *Chemistry::OpenBabelc::OBConversion_ReadString;
*ReadFile = *Chemistry::OpenBabelc::OBConversion_ReadFile;
*BatchFileName = *Chemistry::OpenBabelc::OBConversion_BatchFileName;
*IncrementedFileName = *Chemistry::OpenBabelc::OBConversion_IncrementedFileName;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBResidue ##############

package Chemistry::OpenBabel::OBResidue;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBResidue(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBResidue($self);
        delete $OWNER{$self};
    }
}

*AddAtom = *Chemistry::OpenBabelc::OBResidue_AddAtom;
*InsertAtom = *Chemistry::OpenBabelc::OBResidue_InsertAtom;
*RemoveAtom = *Chemistry::OpenBabelc::OBResidue_RemoveAtom;
*Clear = *Chemistry::OpenBabelc::OBResidue_Clear;
*SetName = *Chemistry::OpenBabelc::OBResidue_SetName;
*SetNum = *Chemistry::OpenBabelc::OBResidue_SetNum;
*SetChain = *Chemistry::OpenBabelc::OBResidue_SetChain;
*SetChainNum = *Chemistry::OpenBabelc::OBResidue_SetChainNum;
*SetIdx = *Chemistry::OpenBabelc::OBResidue_SetIdx;
*SetAtomID = *Chemistry::OpenBabelc::OBResidue_SetAtomID;
*SetHetAtom = *Chemistry::OpenBabelc::OBResidue_SetHetAtom;
*SetSerialNum = *Chemistry::OpenBabelc::OBResidue_SetSerialNum;
*GetName = *Chemistry::OpenBabelc::OBResidue_GetName;
*GetNum = *Chemistry::OpenBabelc::OBResidue_GetNum;
*GetNumAtoms = *Chemistry::OpenBabelc::OBResidue_GetNumAtoms;
*GetChain = *Chemistry::OpenBabelc::OBResidue_GetChain;
*GetChainNum = *Chemistry::OpenBabelc::OBResidue_GetChainNum;
*GetIdx = *Chemistry::OpenBabelc::OBResidue_GetIdx;
*GetResKey = *Chemistry::OpenBabelc::OBResidue_GetResKey;
*GetAtoms = *Chemistry::OpenBabelc::OBResidue_GetAtoms;
*GetBonds = *Chemistry::OpenBabelc::OBResidue_GetBonds;
*GetAtomID = *Chemistry::OpenBabelc::OBResidue_GetAtomID;
*GetSerialNum = *Chemistry::OpenBabelc::OBResidue_GetSerialNum;
*GetAminoAcidProperty = *Chemistry::OpenBabelc::OBResidue_GetAminoAcidProperty;
*GetAtomProperty = *Chemistry::OpenBabelc::OBResidue_GetAtomProperty;
*GetResidueProperty = *Chemistry::OpenBabelc::OBResidue_GetResidueProperty;
*IsHetAtom = *Chemistry::OpenBabelc::OBResidue_IsHetAtom;
*IsResidueType = *Chemistry::OpenBabelc::OBResidue_IsResidueType;
*BeginAtom = *Chemistry::OpenBabelc::OBResidue_BeginAtom;
*NextAtom = *Chemistry::OpenBabelc::OBResidue_NextAtom;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBAtom ##############

package Chemistry::OpenBabel::OBAtom;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBNodeBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBAtom(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBAtom($self);
        delete $OWNER{$self};
    }
}

*Clear = *Chemistry::OpenBabelc::OBAtom_Clear;
*SetIdx = *Chemistry::OpenBabelc::OBAtom_SetIdx;
*SetHyb = *Chemistry::OpenBabelc::OBAtom_SetHyb;
*SetAtomicNum = *Chemistry::OpenBabelc::OBAtom_SetAtomicNum;
*SetIsotope = *Chemistry::OpenBabelc::OBAtom_SetIsotope;
*SetImplicitValence = *Chemistry::OpenBabelc::OBAtom_SetImplicitValence;
*IncrementImplicitValence = *Chemistry::OpenBabelc::OBAtom_IncrementImplicitValence;
*DecrementImplicitValence = *Chemistry::OpenBabelc::OBAtom_DecrementImplicitValence;
*SetFormalCharge = *Chemistry::OpenBabelc::OBAtom_SetFormalCharge;
*SetSpinMultiplicity = *Chemistry::OpenBabelc::OBAtom_SetSpinMultiplicity;
*SetType = *Chemistry::OpenBabelc::OBAtom_SetType;
*SetPartialCharge = *Chemistry::OpenBabelc::OBAtom_SetPartialCharge;
*SetCoordPtr = *Chemistry::OpenBabelc::OBAtom_SetCoordPtr;
*SetVector = *Chemistry::OpenBabelc::OBAtom_SetVector;
*SetResidue = *Chemistry::OpenBabelc::OBAtom_SetResidue;
*SetAromatic = *Chemistry::OpenBabelc::OBAtom_SetAromatic;
*UnsetAromatic = *Chemistry::OpenBabelc::OBAtom_UnsetAromatic;
*SetClockwiseStereo = *Chemistry::OpenBabelc::OBAtom_SetClockwiseStereo;
*SetAntiClockwiseStereo = *Chemistry::OpenBabelc::OBAtom_SetAntiClockwiseStereo;
*SetPositiveStereo = *Chemistry::OpenBabelc::OBAtom_SetPositiveStereo;
*SetNegativeStereo = *Chemistry::OpenBabelc::OBAtom_SetNegativeStereo;
*UnsetStereo = *Chemistry::OpenBabelc::OBAtom_UnsetStereo;
*SetInRing = *Chemistry::OpenBabelc::OBAtom_SetInRing;
*SetChiral = *Chemistry::OpenBabelc::OBAtom_SetChiral;
*ClearCoordPtr = *Chemistry::OpenBabelc::OBAtom_ClearCoordPtr;
*GetFormalCharge = *Chemistry::OpenBabelc::OBAtom_GetFormalCharge;
*GetAtomicNum = *Chemistry::OpenBabelc::OBAtom_GetAtomicNum;
*GetIsotope = *Chemistry::OpenBabelc::OBAtom_GetIsotope;
*GetSpinMultiplicity = *Chemistry::OpenBabelc::OBAtom_GetSpinMultiplicity;
*GetAtomicMass = *Chemistry::OpenBabelc::OBAtom_GetAtomicMass;
*GetExactMass = *Chemistry::OpenBabelc::OBAtom_GetExactMass;
*GetIdx = *Chemistry::OpenBabelc::OBAtom_GetIdx;
*GetCoordinateIdx = *Chemistry::OpenBabelc::OBAtom_GetCoordinateIdx;
*GetCIdx = *Chemistry::OpenBabelc::OBAtom_GetCIdx;
*GetValence = *Chemistry::OpenBabelc::OBAtom_GetValence;
*GetHyb = *Chemistry::OpenBabelc::OBAtom_GetHyb;
*GetImplicitValence = *Chemistry::OpenBabelc::OBAtom_GetImplicitValence;
*GetHvyValence = *Chemistry::OpenBabelc::OBAtom_GetHvyValence;
*GetHeteroValence = *Chemistry::OpenBabelc::OBAtom_GetHeteroValence;
*GetType = *Chemistry::OpenBabelc::OBAtom_GetType;
*GetX = *Chemistry::OpenBabelc::OBAtom_GetX;
*GetY = *Chemistry::OpenBabelc::OBAtom_GetY;
*GetZ = *Chemistry::OpenBabelc::OBAtom_GetZ;
*x = *Chemistry::OpenBabelc::OBAtom_x;
*y = *Chemistry::OpenBabelc::OBAtom_y;
*z = *Chemistry::OpenBabelc::OBAtom_z;
*GetCoordinate = *Chemistry::OpenBabelc::OBAtom_GetCoordinate;
*GetVector = *Chemistry::OpenBabelc::OBAtom_GetVector;
*GetPartialCharge = *Chemistry::OpenBabelc::OBAtom_GetPartialCharge;
*GetResidue = *Chemistry::OpenBabelc::OBAtom_GetResidue;
*GetNewBondVector = *Chemistry::OpenBabelc::OBAtom_GetNewBondVector;
*GetBond = *Chemistry::OpenBabelc::OBAtom_GetBond;
*GetNextAtom = *Chemistry::OpenBabelc::OBAtom_GetNextAtom;
*BeginBonds = *Chemistry::OpenBabelc::OBAtom_BeginBonds;
*EndBonds = *Chemistry::OpenBabelc::OBAtom_EndBonds;
*BeginBond = *Chemistry::OpenBabelc::OBAtom_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBAtom_NextBond;
*BeginNbrAtom = *Chemistry::OpenBabelc::OBAtom_BeginNbrAtom;
*NextNbrAtom = *Chemistry::OpenBabelc::OBAtom_NextNbrAtom;
*GetDistance = *Chemistry::OpenBabelc::OBAtom_GetDistance;
*GetAngle = *Chemistry::OpenBabelc::OBAtom_GetAngle;
*NewResidue = *Chemistry::OpenBabelc::OBAtom_NewResidue;
*DeleteResidue = *Chemistry::OpenBabelc::OBAtom_DeleteResidue;
*AddBond = *Chemistry::OpenBabelc::OBAtom_AddBond;
*InsertBond = *Chemistry::OpenBabelc::OBAtom_InsertBond;
*DeleteBond = *Chemistry::OpenBabelc::OBAtom_DeleteBond;
*ClearBond = *Chemistry::OpenBabelc::OBAtom_ClearBond;
*CountFreeOxygens = *Chemistry::OpenBabelc::OBAtom_CountFreeOxygens;
*ImplicitHydrogenCount = *Chemistry::OpenBabelc::OBAtom_ImplicitHydrogenCount;
*ExplicitHydrogenCount = *Chemistry::OpenBabelc::OBAtom_ExplicitHydrogenCount;
*MemberOfRingCount = *Chemistry::OpenBabelc::OBAtom_MemberOfRingCount;
*MemberOfRingSize = *Chemistry::OpenBabelc::OBAtom_MemberOfRingSize;
*CountRingBonds = *Chemistry::OpenBabelc::OBAtom_CountRingBonds;
*SmallestBondAngle = *Chemistry::OpenBabelc::OBAtom_SmallestBondAngle;
*AverageBondAngle = *Chemistry::OpenBabelc::OBAtom_AverageBondAngle;
*BOSum = *Chemistry::OpenBabelc::OBAtom_BOSum;
*KBOSum = *Chemistry::OpenBabelc::OBAtom_KBOSum;
*HtoMethyl = *Chemistry::OpenBabelc::OBAtom_HtoMethyl;
*SetHybAndGeom = *Chemistry::OpenBabelc::OBAtom_SetHybAndGeom;
*ForceNoH = *Chemistry::OpenBabelc::OBAtom_ForceNoH;
*HasNoHForced = *Chemistry::OpenBabelc::OBAtom_HasNoHForced;
*HasResidue = *Chemistry::OpenBabelc::OBAtom_HasResidue;
*IsHydrogen = *Chemistry::OpenBabelc::OBAtom_IsHydrogen;
*IsCarbon = *Chemistry::OpenBabelc::OBAtom_IsCarbon;
*IsNitrogen = *Chemistry::OpenBabelc::OBAtom_IsNitrogen;
*IsOxygen = *Chemistry::OpenBabelc::OBAtom_IsOxygen;
*IsSulfur = *Chemistry::OpenBabelc::OBAtom_IsSulfur;
*IsPhosphorus = *Chemistry::OpenBabelc::OBAtom_IsPhosphorus;
*IsAromatic = *Chemistry::OpenBabelc::OBAtom_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBAtom_IsInRing;
*IsInRingSize = *Chemistry::OpenBabelc::OBAtom_IsInRingSize;
*IsHeteroatom = *Chemistry::OpenBabelc::OBAtom_IsHeteroatom;
*IsNotCorH = *Chemistry::OpenBabelc::OBAtom_IsNotCorH;
*IsConnected = *Chemistry::OpenBabelc::OBAtom_IsConnected;
*IsOneThree = *Chemistry::OpenBabelc::OBAtom_IsOneThree;
*IsOneFour = *Chemistry::OpenBabelc::OBAtom_IsOneFour;
*IsCarboxylOxygen = *Chemistry::OpenBabelc::OBAtom_IsCarboxylOxygen;
*IsPhosphateOxygen = *Chemistry::OpenBabelc::OBAtom_IsPhosphateOxygen;
*IsSulfateOxygen = *Chemistry::OpenBabelc::OBAtom_IsSulfateOxygen;
*IsNitroOxygen = *Chemistry::OpenBabelc::OBAtom_IsNitroOxygen;
*IsAmideNitrogen = *Chemistry::OpenBabelc::OBAtom_IsAmideNitrogen;
*IsPolarHydrogen = *Chemistry::OpenBabelc::OBAtom_IsPolarHydrogen;
*IsNonPolarHydrogen = *Chemistry::OpenBabelc::OBAtom_IsNonPolarHydrogen;
*IsAromaticNOxide = *Chemistry::OpenBabelc::OBAtom_IsAromaticNOxide;
*IsChiral = *Chemistry::OpenBabelc::OBAtom_IsChiral;
*IsAxial = *Chemistry::OpenBabelc::OBAtom_IsAxial;
*IsClockwise = *Chemistry::OpenBabelc::OBAtom_IsClockwise;
*IsAntiClockwise = *Chemistry::OpenBabelc::OBAtom_IsAntiClockwise;
*IsPositiveStereo = *Chemistry::OpenBabelc::OBAtom_IsPositiveStereo;
*IsNegativeStereo = *Chemistry::OpenBabelc::OBAtom_IsNegativeStereo;
*HasChiralitySpecified = *Chemistry::OpenBabelc::OBAtom_HasChiralitySpecified;
*HasChiralVolume = *Chemistry::OpenBabelc::OBAtom_HasChiralVolume;
*IsHbondAcceptor = *Chemistry::OpenBabelc::OBAtom_IsHbondAcceptor;
*IsHbondDonor = *Chemistry::OpenBabelc::OBAtom_IsHbondDonor;
*IsHbondDonorH = *Chemistry::OpenBabelc::OBAtom_IsHbondDonorH;
*HasAlphaBetaUnsat = *Chemistry::OpenBabelc::OBAtom_HasAlphaBetaUnsat;
*HasBondOfOrder = *Chemistry::OpenBabelc::OBAtom_HasBondOfOrder;
*CountBondsOfOrder = *Chemistry::OpenBabelc::OBAtom_CountBondsOfOrder;
*HasNonSingleBond = *Chemistry::OpenBabelc::OBAtom_HasNonSingleBond;
*HasSingleBond = *Chemistry::OpenBabelc::OBAtom_HasSingleBond;
*HasDoubleBond = *Chemistry::OpenBabelc::OBAtom_HasDoubleBond;
*HasAromaticBond = *Chemistry::OpenBabelc::OBAtom_HasAromaticBond;
*MatchesSMARTS = *Chemistry::OpenBabelc::OBAtom_MatchesSMARTS;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBBond ##############

package Chemistry::OpenBabel::OBBond;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBEdgeBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBBond(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBBond($self);
        delete $OWNER{$self};
    }
}

*SetIdx = *Chemistry::OpenBabelc::OBBond_SetIdx;
*SetBO = *Chemistry::OpenBabelc::OBBond_SetBO;
*SetBegin = *Chemistry::OpenBabelc::OBBond_SetBegin;
*SetEnd = *Chemistry::OpenBabelc::OBBond_SetEnd;
*SetLength = *Chemistry::OpenBabelc::OBBond_SetLength;
*Set = *Chemistry::OpenBabelc::OBBond_Set;
*SetKSingle = *Chemistry::OpenBabelc::OBBond_SetKSingle;
*SetKDouble = *Chemistry::OpenBabelc::OBBond_SetKDouble;
*SetKTriple = *Chemistry::OpenBabelc::OBBond_SetKTriple;
*SetAromatic = *Chemistry::OpenBabelc::OBBond_SetAromatic;
*SetHash = *Chemistry::OpenBabelc::OBBond_SetHash;
*SetWedge = *Chemistry::OpenBabelc::OBBond_SetWedge;
*SetUp = *Chemistry::OpenBabelc::OBBond_SetUp;
*SetDown = *Chemistry::OpenBabelc::OBBond_SetDown;
*SetInRing = *Chemistry::OpenBabelc::OBBond_SetInRing;
*SetClosure = *Chemistry::OpenBabelc::OBBond_SetClosure;
*UnsetHash = *Chemistry::OpenBabelc::OBBond_UnsetHash;
*UnsetWedge = *Chemistry::OpenBabelc::OBBond_UnsetWedge;
*UnsetUp = *Chemistry::OpenBabelc::OBBond_UnsetUp;
*UnsetDown = *Chemistry::OpenBabelc::OBBond_UnsetDown;
*UnsetAromatic = *Chemistry::OpenBabelc::OBBond_UnsetAromatic;
*UnsetKekule = *Chemistry::OpenBabelc::OBBond_UnsetKekule;
*GetBO = *Chemistry::OpenBabelc::OBBond_GetBO;
*GetBondOrder = *Chemistry::OpenBabelc::OBBond_GetBondOrder;
*GetFlags = *Chemistry::OpenBabelc::OBBond_GetFlags;
*GetBeginAtomIdx = *Chemistry::OpenBabelc::OBBond_GetBeginAtomIdx;
*GetEndAtomIdx = *Chemistry::OpenBabelc::OBBond_GetEndAtomIdx;
*GetBeginAtom = *Chemistry::OpenBabelc::OBBond_GetBeginAtom;
*GetEndAtom = *Chemistry::OpenBabelc::OBBond_GetEndAtom;
*GetNbrAtom = *Chemistry::OpenBabelc::OBBond_GetNbrAtom;
*GetEquibLength = *Chemistry::OpenBabelc::OBBond_GetEquibLength;
*GetLength = *Chemistry::OpenBabelc::OBBond_GetLength;
*GetNbrAtomIdx = *Chemistry::OpenBabelc::OBBond_GetNbrAtomIdx;
*IsAromatic = *Chemistry::OpenBabelc::OBBond_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBBond_IsInRing;
*IsRotor = *Chemistry::OpenBabelc::OBBond_IsRotor;
*IsAmide = *Chemistry::OpenBabelc::OBBond_IsAmide;
*IsPrimaryAmide = *Chemistry::OpenBabelc::OBBond_IsPrimaryAmide;
*IsSecondaryAmide = *Chemistry::OpenBabelc::OBBond_IsSecondaryAmide;
*IsEster = *Chemistry::OpenBabelc::OBBond_IsEster;
*IsCarbonyl = *Chemistry::OpenBabelc::OBBond_IsCarbonyl;
*IsSingle = *Chemistry::OpenBabelc::OBBond_IsSingle;
*IsDouble = *Chemistry::OpenBabelc::OBBond_IsDouble;
*IsTriple = *Chemistry::OpenBabelc::OBBond_IsTriple;
*IsKSingle = *Chemistry::OpenBabelc::OBBond_IsKSingle;
*IsKDouble = *Chemistry::OpenBabelc::OBBond_IsKDouble;
*IsKTriple = *Chemistry::OpenBabelc::OBBond_IsKTriple;
*IsClosure = *Chemistry::OpenBabelc::OBBond_IsClosure;
*IsUp = *Chemistry::OpenBabelc::OBBond_IsUp;
*IsDown = *Chemistry::OpenBabelc::OBBond_IsDown;
*IsWedge = *Chemistry::OpenBabelc::OBBond_IsWedge;
*IsHash = *Chemistry::OpenBabelc::OBBond_IsHash;
*IsDoubleBondGeometry = *Chemistry::OpenBabelc::OBBond_IsDoubleBondGeometry;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBMol ##############

package Chemistry::OpenBabel::OBMol;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel::OBGraphBase Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBMol(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBMol($self);
        delete $OWNER{$self};
    }
}

*ReserveAtoms = *Chemistry::OpenBabelc::OBMol_ReserveAtoms;
*CreateAtom = *Chemistry::OpenBabelc::OBMol_CreateAtom;
*CreateBond = *Chemistry::OpenBabelc::OBMol_CreateBond;
*DestroyAtom = *Chemistry::OpenBabelc::OBMol_DestroyAtom;
*DestroyBond = *Chemistry::OpenBabelc::OBMol_DestroyBond;
*AddAtom = *Chemistry::OpenBabelc::OBMol_AddAtom;
*AddBond = *Chemistry::OpenBabelc::OBMol_AddBond;
*AddResidue = *Chemistry::OpenBabelc::OBMol_AddResidue;
*InsertAtom = *Chemistry::OpenBabelc::OBMol_InsertAtom;
*DeleteAtom = *Chemistry::OpenBabelc::OBMol_DeleteAtom;
*DeleteBond = *Chemistry::OpenBabelc::OBMol_DeleteBond;
*DeleteResidue = *Chemistry::OpenBabelc::OBMol_DeleteResidue;
*NewAtom = *Chemistry::OpenBabelc::OBMol_NewAtom;
*NewResidue = *Chemistry::OpenBabelc::OBMol_NewResidue;
*BeginModify = *Chemistry::OpenBabelc::OBMol_BeginModify;
*EndModify = *Chemistry::OpenBabelc::OBMol_EndModify;
*GetMod = *Chemistry::OpenBabelc::OBMol_GetMod;
*IncrementMod = *Chemistry::OpenBabelc::OBMol_IncrementMod;
*DecrementMod = *Chemistry::OpenBabelc::OBMol_DecrementMod;
*GetFlags = *Chemistry::OpenBabelc::OBMol_GetFlags;
*GetTitle = *Chemistry::OpenBabelc::OBMol_GetTitle;
*NumAtoms = *Chemistry::OpenBabelc::OBMol_NumAtoms;
*NumBonds = *Chemistry::OpenBabelc::OBMol_NumBonds;
*NumHvyAtoms = *Chemistry::OpenBabelc::OBMol_NumHvyAtoms;
*NumResidues = *Chemistry::OpenBabelc::OBMol_NumResidues;
*NumRotors = *Chemistry::OpenBabelc::OBMol_NumRotors;
*GetAtom = *Chemistry::OpenBabelc::OBMol_GetAtom;
*GetFirstAtom = *Chemistry::OpenBabelc::OBMol_GetFirstAtom;
*GetBond = *Chemistry::OpenBabelc::OBMol_GetBond;
*GetResidue = *Chemistry::OpenBabelc::OBMol_GetResidue;
*GetInternalCoord = *Chemistry::OpenBabelc::OBMol_GetInternalCoord;
*GetTorsion = *Chemistry::OpenBabelc::OBMol_GetTorsion;
*GetFormula = *Chemistry::OpenBabelc::OBMol_GetFormula;
*GetSpacedFormula = *Chemistry::OpenBabelc::OBMol_GetSpacedFormula;
*GetEnergy = *Chemistry::OpenBabelc::OBMol_GetEnergy;
*GetMolWt = *Chemistry::OpenBabelc::OBMol_GetMolWt;
*GetExactMass = *Chemistry::OpenBabelc::OBMol_GetExactMass;
*GetTotalCharge = *Chemistry::OpenBabelc::OBMol_GetTotalCharge;
*GetTotalSpinMultiplicity = *Chemistry::OpenBabelc::OBMol_GetTotalSpinMultiplicity;
*GetDimension = *Chemistry::OpenBabelc::OBMol_GetDimension;
*GetCoordinates = *Chemistry::OpenBabelc::OBMol_GetCoordinates;
*GetSSSR = *Chemistry::OpenBabelc::OBMol_GetSSSR;
*AutomaticFormalCharge = *Chemistry::OpenBabelc::OBMol_AutomaticFormalCharge;
*AutomaticPartialCharge = *Chemistry::OpenBabelc::OBMol_AutomaticPartialCharge;
*SetTitle = *Chemistry::OpenBabelc::OBMol_SetTitle;
*SetFormula = *Chemistry::OpenBabelc::OBMol_SetFormula;
*SetEnergy = *Chemistry::OpenBabelc::OBMol_SetEnergy;
*SetDimension = *Chemistry::OpenBabelc::OBMol_SetDimension;
*SetTotalCharge = *Chemistry::OpenBabelc::OBMol_SetTotalCharge;
*SetTotalSpinMultiplicity = *Chemistry::OpenBabelc::OBMol_SetTotalSpinMultiplicity;
*SetInternalCoord = *Chemistry::OpenBabelc::OBMol_SetInternalCoord;
*SetAutomaticFormalCharge = *Chemistry::OpenBabelc::OBMol_SetAutomaticFormalCharge;
*SetAutomaticPartialCharge = *Chemistry::OpenBabelc::OBMol_SetAutomaticPartialCharge;
*SetAromaticPerceived = *Chemistry::OpenBabelc::OBMol_SetAromaticPerceived;
*SetSSSRPerceived = *Chemistry::OpenBabelc::OBMol_SetSSSRPerceived;
*SetRingAtomsAndBondsPerceived = *Chemistry::OpenBabelc::OBMol_SetRingAtomsAndBondsPerceived;
*SetAtomTypesPerceived = *Chemistry::OpenBabelc::OBMol_SetAtomTypesPerceived;
*SetChainsPerceived = *Chemistry::OpenBabelc::OBMol_SetChainsPerceived;
*SetChiralityPerceived = *Chemistry::OpenBabelc::OBMol_SetChiralityPerceived;
*SetPartialChargesPerceived = *Chemistry::OpenBabelc::OBMol_SetPartialChargesPerceived;
*SetHybridizationPerceived = *Chemistry::OpenBabelc::OBMol_SetHybridizationPerceived;
*SetImplicitValencePerceived = *Chemistry::OpenBabelc::OBMol_SetImplicitValencePerceived;
*SetKekulePerceived = *Chemistry::OpenBabelc::OBMol_SetKekulePerceived;
*SetClosureBondsPerceived = *Chemistry::OpenBabelc::OBMol_SetClosureBondsPerceived;
*SetHydrogensAdded = *Chemistry::OpenBabelc::OBMol_SetHydrogensAdded;
*SetCorrectedForPH = *Chemistry::OpenBabelc::OBMol_SetCorrectedForPH;
*SetAromaticCorrected = *Chemistry::OpenBabelc::OBMol_SetAromaticCorrected;
*SetSpinMultiplicityAssigned = *Chemistry::OpenBabelc::OBMol_SetSpinMultiplicityAssigned;
*SetFlags = *Chemistry::OpenBabelc::OBMol_SetFlags;
*UnsetAromaticPerceived = *Chemistry::OpenBabelc::OBMol_UnsetAromaticPerceived;
*UnsetPartialChargesPerceived = *Chemistry::OpenBabelc::OBMol_UnsetPartialChargesPerceived;
*UnsetImplicitValencePerceived = *Chemistry::OpenBabelc::OBMol_UnsetImplicitValencePerceived;
*UnsetFlag = *Chemistry::OpenBabelc::OBMol_UnsetFlag;
*DoTransformations = *Chemistry::OpenBabelc::OBMol_DoTransformations;
*ClassDescription = *Chemistry::OpenBabelc::OBMol_ClassDescription;
*Clear = *Chemistry::OpenBabelc::OBMol_Clear;
*RenumberAtoms = *Chemistry::OpenBabelc::OBMol_RenumberAtoms;
*ToInertialFrame = *Chemistry::OpenBabelc::OBMol_ToInertialFrame;
*Translate = *Chemistry::OpenBabelc::OBMol_Translate;
*Rotate = *Chemistry::OpenBabelc::OBMol_Rotate;
*Kekulize = *Chemistry::OpenBabelc::OBMol_Kekulize;
*PerceiveKekuleBonds = *Chemistry::OpenBabelc::OBMol_PerceiveKekuleBonds;
*NewPerceiveKekuleBonds = *Chemistry::OpenBabelc::OBMol_NewPerceiveKekuleBonds;
*DeleteHydrogen = *Chemistry::OpenBabelc::OBMol_DeleteHydrogen;
*DeleteHydrogens = *Chemistry::OpenBabelc::OBMol_DeleteHydrogens;
*DeleteNonPolarHydrogens = *Chemistry::OpenBabelc::OBMol_DeleteNonPolarHydrogens;
*AddHydrogens = *Chemistry::OpenBabelc::OBMol_AddHydrogens;
*AddPolarHydrogens = *Chemistry::OpenBabelc::OBMol_AddPolarHydrogens;
*StripSalts = *Chemistry::OpenBabelc::OBMol_StripSalts;
*ConvertDativeBonds = *Chemistry::OpenBabelc::OBMol_ConvertDativeBonds;
*CorrectForPH = *Chemistry::OpenBabelc::OBMol_CorrectForPH;
*AssignSpinMultiplicity = *Chemistry::OpenBabelc::OBMol_AssignSpinMultiplicity;
*Center = *Chemistry::OpenBabelc::OBMol_Center;
*SetTorsion = *Chemistry::OpenBabelc::OBMol_SetTorsion;
*FindSSSR = *Chemistry::OpenBabelc::OBMol_FindSSSR;
*FindRingAtomsAndBonds = *Chemistry::OpenBabelc::OBMol_FindRingAtomsAndBonds;
*FindChiralCenters = *Chemistry::OpenBabelc::OBMol_FindChiralCenters;
*FindChildren = *Chemistry::OpenBabelc::OBMol_FindChildren;
*FindLargestFragment = *Chemistry::OpenBabelc::OBMol_FindLargestFragment;
*ContigFragList = *Chemistry::OpenBabelc::OBMol_ContigFragList;
*Align = *Chemistry::OpenBabelc::OBMol_Align;
*ConnectTheDots = *Chemistry::OpenBabelc::OBMol_ConnectTheDots;
*PerceiveBondOrders = *Chemistry::OpenBabelc::OBMol_PerceiveBondOrders;
*FindTorsions = *Chemistry::OpenBabelc::OBMol_FindTorsions;
*GetGTDVector = *Chemistry::OpenBabelc::OBMol_GetGTDVector;
*GetGIVector = *Chemistry::OpenBabelc::OBMol_GetGIVector;
*GetGIDVector = *Chemistry::OpenBabelc::OBMol_GetGIDVector;
*Has2D = *Chemistry::OpenBabelc::OBMol_Has2D;
*Has3D = *Chemistry::OpenBabelc::OBMol_Has3D;
*HasNonZeroCoords = *Chemistry::OpenBabelc::OBMol_HasNonZeroCoords;
*HasAromaticPerceived = *Chemistry::OpenBabelc::OBMol_HasAromaticPerceived;
*HasSSSRPerceived = *Chemistry::OpenBabelc::OBMol_HasSSSRPerceived;
*HasRingAtomsAndBondsPerceived = *Chemistry::OpenBabelc::OBMol_HasRingAtomsAndBondsPerceived;
*HasAtomTypesPerceived = *Chemistry::OpenBabelc::OBMol_HasAtomTypesPerceived;
*HasChiralityPerceived = *Chemistry::OpenBabelc::OBMol_HasChiralityPerceived;
*HasPartialChargesPerceived = *Chemistry::OpenBabelc::OBMol_HasPartialChargesPerceived;
*HasHybridizationPerceived = *Chemistry::OpenBabelc::OBMol_HasHybridizationPerceived;
*HasImplicitValencePerceived = *Chemistry::OpenBabelc::OBMol_HasImplicitValencePerceived;
*HasKekulePerceived = *Chemistry::OpenBabelc::OBMol_HasKekulePerceived;
*HasClosureBondsPerceived = *Chemistry::OpenBabelc::OBMol_HasClosureBondsPerceived;
*HasChainsPerceived = *Chemistry::OpenBabelc::OBMol_HasChainsPerceived;
*HasHydrogensAdded = *Chemistry::OpenBabelc::OBMol_HasHydrogensAdded;
*HasAromaticCorrected = *Chemistry::OpenBabelc::OBMol_HasAromaticCorrected;
*IsCorrectedForPH = *Chemistry::OpenBabelc::OBMol_IsCorrectedForPH;
*HasSpinMultiplicityAssigned = *Chemistry::OpenBabelc::OBMol_HasSpinMultiplicityAssigned;
*IsChiral = *Chemistry::OpenBabelc::OBMol_IsChiral;
*Empty = *Chemistry::OpenBabelc::OBMol_Empty;
*NumConformers = *Chemistry::OpenBabelc::OBMol_NumConformers;
*SetConformers = *Chemistry::OpenBabelc::OBMol_SetConformers;
*AddConformer = *Chemistry::OpenBabelc::OBMol_AddConformer;
*SetConformer = *Chemistry::OpenBabelc::OBMol_SetConformer;
*CopyConformer = *Chemistry::OpenBabelc::OBMol_CopyConformer;
*DeleteConformer = *Chemistry::OpenBabelc::OBMol_DeleteConformer;
*GetConformer = *Chemistry::OpenBabelc::OBMol_GetConformer;
*BeginConformer = *Chemistry::OpenBabelc::OBMol_BeginConformer;
*NextConformer = *Chemistry::OpenBabelc::OBMol_NextConformer;
*GetConformers = *Chemistry::OpenBabelc::OBMol_GetConformers;
*BeginAtom = *Chemistry::OpenBabelc::OBMol_BeginAtom;
*NextAtom = *Chemistry::OpenBabelc::OBMol_NextAtom;
*BeginBond = *Chemistry::OpenBabelc::OBMol_BeginBond;
*NextBond = *Chemistry::OpenBabelc::OBMol_NextBond;
*BeginResidue = *Chemistry::OpenBabelc::OBMol_BeginResidue;
*NextResidue = *Chemistry::OpenBabelc::OBMol_NextResidue;
*BeginInternalCoord = *Chemistry::OpenBabelc::OBMol_BeginInternalCoord;
*NextInternalCoord = *Chemistry::OpenBabelc::OBMol_NextInternalCoord;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBInternalCoord ##############

package Chemistry::OpenBabel::OBInternalCoord;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig__a_get = *Chemistry::OpenBabelc::OBInternalCoord__a_get;
*swig__a_set = *Chemistry::OpenBabelc::OBInternalCoord__a_set;
*swig__b_get = *Chemistry::OpenBabelc::OBInternalCoord__b_get;
*swig__b_set = *Chemistry::OpenBabelc::OBInternalCoord__b_set;
*swig__c_get = *Chemistry::OpenBabelc::OBInternalCoord__c_get;
*swig__c_set = *Chemistry::OpenBabelc::OBInternalCoord__c_set;
*swig__dst_get = *Chemistry::OpenBabelc::OBInternalCoord__dst_get;
*swig__dst_set = *Chemistry::OpenBabelc::OBInternalCoord__dst_set;
*swig__ang_get = *Chemistry::OpenBabelc::OBInternalCoord__ang_get;
*swig__ang_set = *Chemistry::OpenBabelc::OBInternalCoord__ang_set;
*swig__tor_get = *Chemistry::OpenBabelc::OBInternalCoord__tor_get;
*swig__tor_set = *Chemistry::OpenBabelc::OBInternalCoord__tor_set;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBInternalCoord(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBInternalCoord($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRTree ##############

package Chemistry::OpenBabel::OBRTree;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRTree(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRTree($self);
        delete $OWNER{$self};
    }
}

*GetAtomIdx = *Chemistry::OpenBabelc::OBRTree_GetAtomIdx;
*PathToRoot = *Chemistry::OpenBabelc::OBRTree_PathToRoot;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRing ##############

package Chemistry::OpenBabel::OBRing;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
*swig__path_get = *Chemistry::OpenBabelc::OBRing__path_get;
*swig__path_set = *Chemistry::OpenBabelc::OBRing__path_set;
*swig__pathset_get = *Chemistry::OpenBabelc::OBRing__pathset_get;
*swig__pathset_set = *Chemistry::OpenBabelc::OBRing__pathset_set;
*findCenterAndNormal = *Chemistry::OpenBabelc::OBRing_findCenterAndNormal;
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRing(@_);
    bless $self, $pkg if defined($self);
}

*Size = *Chemistry::OpenBabelc::OBRing_Size;
*PathSize = *Chemistry::OpenBabelc::OBRing_PathSize;
*IsMember = *Chemistry::OpenBabelc::OBRing_IsMember;
*IsAromatic = *Chemistry::OpenBabelc::OBRing_IsAromatic;
*IsInRing = *Chemistry::OpenBabelc::OBRing_IsInRing;
*SetParent = *Chemistry::OpenBabelc::OBRing_SetParent;
*GetParent = *Chemistry::OpenBabelc::OBRing_GetParent;
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRing($self);
        delete $OWNER{$self};
    }
}

sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBRingSearch ##############

package Chemistry::OpenBabel::OBRingSearch;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBRingSearch(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBRingSearch($self);
        delete $OWNER{$self};
    }
}

*SortRings = *Chemistry::OpenBabelc::OBRingSearch_SortRings;
*RemoveRedundant = *Chemistry::OpenBabelc::OBRingSearch_RemoveRedundant;
*AddRingFromClosure = *Chemistry::OpenBabelc::OBRingSearch_AddRingFromClosure;
*WriteRings = *Chemistry::OpenBabelc::OBRingSearch_WriteRings;
*SaveUniqueRing = *Chemistry::OpenBabelc::OBRingSearch_SaveUniqueRing;
*BeginRings = *Chemistry::OpenBabelc::OBRingSearch_BeginRings;
*EndRings = *Chemistry::OpenBabelc::OBRingSearch_EndRings;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSmartsPattern ##############

package Chemistry::OpenBabel::OBSmartsPattern;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSmartsPattern($self);
        delete $OWNER{$self};
    }
}

sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSmartsPattern(@_);
    bless $self, $pkg if defined($self);
}

*NumMatches = *Chemistry::OpenBabelc::OBSmartsPattern_NumMatches;
*NumAtoms = *Chemistry::OpenBabelc::OBSmartsPattern_NumAtoms;
*NumBonds = *Chemistry::OpenBabelc::OBSmartsPattern_NumBonds;
*GetAtomicNum = *Chemistry::OpenBabelc::OBSmartsPattern_GetAtomicNum;
*GetBond = *Chemistry::OpenBabelc::OBSmartsPattern_GetBond;
*GetCharge = *Chemistry::OpenBabelc::OBSmartsPattern_GetCharge;
*GetSMARTS = *Chemistry::OpenBabelc::OBSmartsPattern_GetSMARTS;
*GetVectorBinding = *Chemistry::OpenBabelc::OBSmartsPattern_GetVectorBinding;
*Empty = *Chemistry::OpenBabelc::OBSmartsPattern_Empty;
*IsValid = *Chemistry::OpenBabelc::OBSmartsPattern_IsValid;
*Init = *Chemistry::OpenBabelc::OBSmartsPattern_Init;
*WriteMapList = *Chemistry::OpenBabelc::OBSmartsPattern_WriteMapList;
*Match = *Chemistry::OpenBabelc::OBSmartsPattern_Match;
*RestrictedMatch = *Chemistry::OpenBabelc::OBSmartsPattern_RestrictedMatch;
*GetMapList = *Chemistry::OpenBabelc::OBSmartsPattern_GetMapList;
*GetUMapList = *Chemistry::OpenBabelc::OBSmartsPattern_GetUMapList;
*BeginMList = *Chemistry::OpenBabelc::OBSmartsPattern_BeginMList;
*EndMList = *Chemistry::OpenBabelc::OBSmartsPattern_EndMList;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


############# Class : Chemistry::OpenBabel::OBSSMatch ##############

package Chemistry::OpenBabel::OBSSMatch;
use vars qw(@ISA %OWNER %ITERATORS %BLESSEDMEMBERS);
@ISA = qw( Chemistry::OpenBabel );
%OWNER = ();
%ITERATORS = ();
sub new {
    my $pkg = shift;
    my $self = Chemistry::OpenBabelc::new_OBSSMatch(@_);
    bless $self, $pkg if defined($self);
}

sub DESTROY {
    return unless $_[0]->isa('HASH');
    my $self = tied(%{$_[0]});
    return unless defined $self;
    delete $ITERATORS{$self};
    if (exists $OWNER{$self}) {
        Chemistry::OpenBabelc::delete_OBSSMatch($self);
        delete $OWNER{$self};
    }
}

*Match = *Chemistry::OpenBabelc::OBSSMatch_Match;
sub DISOWN {
    my $self = shift;
    my $ptr = tied(%$self);
    delete $OWNER{$ptr};
}

sub ACQUIRE {
    my $self = shift;
    my $ptr = tied(%$self);
    $OWNER{$ptr} = 1;
}


# ------- VARIABLE STUBS --------

package Chemistry::OpenBabel;

*FILE_SEP_CHAR = *Chemistry::OpenBabelc::FILE_SEP_CHAR;
*PI = *Chemistry::OpenBabelc::PI;
*RAD_TO_DEG = *Chemistry::OpenBabelc::RAD_TO_DEG;
*DEG_TO_RAD = *Chemistry::OpenBabelc::DEG_TO_RAD;

my %__VZero_hash;
tie %__VZero_hash,"Chemistry::OpenBabel::vector3", $Chemistry::OpenBabelc::VZero;
$VZero= \%__VZero_hash;
bless $VZero, Chemistry::OpenBabel::vector3;

my %__VX_hash;
tie %__VX_hash,"Chemistry::OpenBabel::vector3", $Chemistry::OpenBabelc::VX;
$VX= \%__VX_hash;
bless $VX, Chemistry::OpenBabel::vector3;

my %__VY_hash;
tie %__VY_hash,"Chemistry::OpenBabel::vector3", $Chemistry::OpenBabelc::VY;
$VY= \%__VY_hash;
bless $VY, Chemistry::OpenBabel::vector3;

my %__VZ_hash;
tie %__VZ_hash,"Chemistry::OpenBabel::vector3", $Chemistry::OpenBabelc::VZ;
$VZ= \%__VZ_hash;
bless $VZ, Chemistry::OpenBabel::vector3;
*UndefinedData = *Chemistry::OpenBabelc::UndefinedData;
*PairData = *Chemistry::OpenBabelc::PairData;
*EnergyData = *Chemistry::OpenBabelc::EnergyData;
*CommentData = *Chemistry::OpenBabelc::CommentData;
*ConformerData = *Chemistry::OpenBabelc::ConformerData;
*ExternalBondData = *Chemistry::OpenBabelc::ExternalBondData;
*RotamerList = *Chemistry::OpenBabelc::RotamerList;
*VirtualBondData = *Chemistry::OpenBabelc::VirtualBondData;
*RingData = *Chemistry::OpenBabelc::RingData;
*TorsionData = *Chemistry::OpenBabelc::TorsionData;
*AngleData = *Chemistry::OpenBabelc::AngleData;
*SerialNums = *Chemistry::OpenBabelc::SerialNums;
*UnitCell = *Chemistry::OpenBabelc::UnitCell;
*SpinData = *Chemistry::OpenBabelc::SpinData;
*ChargeData = *Chemistry::OpenBabelc::ChargeData;
*SymmetryData = *Chemistry::OpenBabelc::SymmetryData;
*ChiralData = *Chemistry::OpenBabelc::ChiralData;
*OccupationData = *Chemistry::OpenBabelc::OccupationData;
*DensityData = *Chemistry::OpenBabelc::DensityData;
*ElectronicData = *Chemistry::OpenBabelc::ElectronicData;
*VibrationData = *Chemistry::OpenBabelc::VibrationData;
*RotationData = *Chemistry::OpenBabelc::RotationData;
*NuclearData = *Chemistry::OpenBabelc::NuclearData;
*CustomData0 = *Chemistry::OpenBabelc::CustomData0;
*CustomData1 = *Chemistry::OpenBabelc::CustomData1;
*CustomData2 = *Chemistry::OpenBabelc::CustomData2;
*CustomData3 = *Chemistry::OpenBabelc::CustomData3;
*CustomData4 = *Chemistry::OpenBabelc::CustomData4;
*CustomData5 = *Chemistry::OpenBabelc::CustomData5;
*CustomData6 = *Chemistry::OpenBabelc::CustomData6;
*CustomData7 = *Chemistry::OpenBabelc::CustomData7;
*CustomData8 = *Chemistry::OpenBabelc::CustomData8;
*CustomData9 = *Chemistry::OpenBabelc::CustomData9;
*CustomData10 = *Chemistry::OpenBabelc::CustomData10;
*CustomData11 = *Chemistry::OpenBabelc::CustomData11;
*CustomData12 = *Chemistry::OpenBabelc::CustomData12;
*CustomData13 = *Chemistry::OpenBabelc::CustomData13;
*CustomData14 = *Chemistry::OpenBabelc::CustomData14;
*CustomData15 = *Chemistry::OpenBabelc::CustomData15;
*output = *Chemistry::OpenBabelc::output;
*input = *Chemistry::OpenBabelc::input;
*calcvolume = *Chemistry::OpenBabelc::calcvolume;
*NOTREADABLE = *Chemistry::OpenBabelc::NOTREADABLE;
*READONEONLY = *Chemistry::OpenBabelc::READONEONLY;
*READBINARY = *Chemistry::OpenBabelc::READBINARY;
*ZEROATOMSOK = *Chemistry::OpenBabelc::ZEROATOMSOK;
*NOTWRITABLE = *Chemistry::OpenBabelc::NOTWRITABLE;
*WRITEONEONLY = *Chemistry::OpenBabelc::WRITEONEONLY;
*WRITEBINARY = *Chemistry::OpenBabelc::WRITEBINARY;
*DEFAULTFORMAT = *Chemistry::OpenBabelc::DEFAULTFORMAT;
*OB_4RING_ATOM = *Chemistry::OpenBabelc::OB_4RING_ATOM;
*OB_3RING_ATOM = *Chemistry::OpenBabelc::OB_3RING_ATOM;
*OB_AROMATIC_ATOM = *Chemistry::OpenBabelc::OB_AROMATIC_ATOM;
*OB_RING_ATOM = *Chemistry::OpenBabelc::OB_RING_ATOM;
*OB_CSTEREO_ATOM = *Chemistry::OpenBabelc::OB_CSTEREO_ATOM;
*OB_ACSTEREO_ATOM = *Chemistry::OpenBabelc::OB_ACSTEREO_ATOM;
*OB_DONOR_ATOM = *Chemistry::OpenBabelc::OB_DONOR_ATOM;
*OB_ACCEPTOR_ATOM = *Chemistry::OpenBabelc::OB_ACCEPTOR_ATOM;
*OB_CHIRAL_ATOM = *Chemistry::OpenBabelc::OB_CHIRAL_ATOM;
*OB_POS_CHIRAL_ATOM = *Chemistry::OpenBabelc::OB_POS_CHIRAL_ATOM;
*OB_NEG_CHIRAL_ATOM = *Chemistry::OpenBabelc::OB_NEG_CHIRAL_ATOM;
*OB_ATOM_HAS_NO_H = *Chemistry::OpenBabelc::OB_ATOM_HAS_NO_H;
*OB_AROMATIC_BOND = *Chemistry::OpenBabelc::OB_AROMATIC_BOND;
*OB_WEDGE_BOND = *Chemistry::OpenBabelc::OB_WEDGE_BOND;
*OB_HASH_BOND = *Chemistry::OpenBabelc::OB_HASH_BOND;
*OB_RING_BOND = *Chemistry::OpenBabelc::OB_RING_BOND;
*OB_TORUP_BOND = *Chemistry::OpenBabelc::OB_TORUP_BOND;
*OB_TORDOWN_BOND = *Chemistry::OpenBabelc::OB_TORDOWN_BOND;
*OB_KSINGLE_BOND = *Chemistry::OpenBabelc::OB_KSINGLE_BOND;
*OB_KDOUBLE_BOND = *Chemistry::OpenBabelc::OB_KDOUBLE_BOND;
*OB_KTRIPLE_BOND = *Chemistry::OpenBabelc::OB_KTRIPLE_BOND;
*OB_CLOSURE_BOND = *Chemistry::OpenBabelc::OB_CLOSURE_BOND;
*OB_SSSR_MOL = *Chemistry::OpenBabelc::OB_SSSR_MOL;
*OB_RINGFLAGS_MOL = *Chemistry::OpenBabelc::OB_RINGFLAGS_MOL;
*OB_AROMATIC_MOL = *Chemistry::OpenBabelc::OB_AROMATIC_MOL;
*OB_ATOMTYPES_MOL = *Chemistry::OpenBabelc::OB_ATOMTYPES_MOL;
*OB_CHIRALITY_MOL = *Chemistry::OpenBabelc::OB_CHIRALITY_MOL;
*OB_PCHARGE_MOL = *Chemistry::OpenBabelc::OB_PCHARGE_MOL;
*OB_HYBRID_MOL = *Chemistry::OpenBabelc::OB_HYBRID_MOL;
*OB_IMPVAL_MOL = *Chemistry::OpenBabelc::OB_IMPVAL_MOL;
*OB_KEKULE_MOL = *Chemistry::OpenBabelc::OB_KEKULE_MOL;
*OB_CLOSURE_MOL = *Chemistry::OpenBabelc::OB_CLOSURE_MOL;
*OB_H_ADDED_MOL = *Chemistry::OpenBabelc::OB_H_ADDED_MOL;
*OB_PH_CORRECTED_MOL = *Chemistry::OpenBabelc::OB_PH_CORRECTED_MOL;
*OB_AROM_CORRECTED_MOL = *Chemistry::OpenBabelc::OB_AROM_CORRECTED_MOL;
*OB_CHAINS_MOL = *Chemistry::OpenBabelc::OB_CHAINS_MOL;
*OB_TCHARGE_MOL = *Chemistry::OpenBabelc::OB_TCHARGE_MOL;
*OB_TSPIN_MOL = *Chemistry::OpenBabelc::OB_TSPIN_MOL;
*OB_CURRENT_CONFORMER = *Chemistry::OpenBabelc::OB_CURRENT_CONFORMER;

my %__etab_hash;
tie %__etab_hash,"Chemistry::OpenBabel::OBElementTable", $Chemistry::OpenBabelc::etab;
$etab= \%__etab_hash;
bless $etab, Chemistry::OpenBabel::OBElementTable;

my %__ttab_hash;
tie %__ttab_hash,"Chemistry::OpenBabel::OBTypeTable", $Chemistry::OpenBabelc::ttab;
$ttab= \%__ttab_hash;
bless $ttab, Chemistry::OpenBabel::OBTypeTable;

my %__isotab_hash;
tie %__isotab_hash,"Chemistry::OpenBabel::OBIsotopeTable", $Chemistry::OpenBabelc::isotab;
$isotab= \%__isotab_hash;
bless $isotab, Chemistry::OpenBabel::OBIsotopeTable;
*aromtyper = *Chemistry::OpenBabelc::aromtyper;
*atomtyper = *Chemistry::OpenBabelc::atomtyper;
*chainsparser = *Chemistry::OpenBabelc::chainsparser;

my %__resdat_hash;
tie %__resdat_hash,"Chemistry::OpenBabel::OBResidueData", $Chemistry::OpenBabelc::resdat;
$resdat= \%__resdat_hash;
bless $resdat, Chemistry::OpenBabel::OBResidueData;
*BUFF_SIZE = *Chemistry::OpenBabelc::BUFF_SIZE;
1;
