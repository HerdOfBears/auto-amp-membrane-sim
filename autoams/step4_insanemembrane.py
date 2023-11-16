# the class needs to take all the possible options for the insane script, then use the os package to run the code
import os
import subprocess

__all__ = ["Insanity"]


class Insanity():

    def __init__(self, sys, Defaults=True):
        self.add_args = None
        self.altail_ = None
        self.allink_ = None
        self.alhead = None
        self.alname_ = None
        self.charge_ = None
        self.salt_ = None
        self.excl_ = None
        self.solr_ = None
        self.sold_ = None
        self.sol_ = None
        self.prot_args = None
        self.dm_ = None
        self.ring_ = None
        self.fudge_ = None
        self.op_ = None
        self.od_ = None
        self.rotate_ = None
        self.orient_ = None
        self.center_ = None
        self.bd_ = None
        self.rand_ = None
        self.disc_ = None
        self.hole_ = None
        self.asym_ = None
        self.au_ = None
        self.a_ = None
        self.u_ = None
        self.l_ = None
        self.mem_args = None
        self.pbc_args = None
        self.pbc_ = None
        self.d_ = None
        self.xyz = None
        self.box_ = None
        self.command = ''
        self.PBS = False
        self.MS = False
        self.PS = False
        self.SlS = False
        self.AO = False
        self.system = sys
        # check for insane.py in the folder
        if os.path.isfile('insane.py'):
            pass
        else:
            raise Exception("Missing script file: insane.py. Please locate and insert into directory")

        print('Initialized insane script object. If you need help, Please run the therapy method!')

    def PeriodicBoundarySetup(self, pbc=None, d=None, dz=None, x=None, y=None, z=None, box=None, n=None):
        """collect kwargs and format into command"""

        # PERIODIC BOUNDARY CONDITIONS---------------------------
        periodic_options = ['hexagonal', 'rectangular', 'square', 'cubic', 'optimal', 'keep']
        if pbc in periodic_options:
            self.pbc_ = "-pbc " + pbc + " "
        else:
            raise Exception(
                "pbc boolean MUST be type 'string' and MUST match one of the accepted options: \nhexagonal, rectangular, square, cubic, optimal or keep")

        # DISTANCE BETWEEN PERIODIC IMAGES ----------------------
        if d == None:
            self.d_ = ""  # dont add the -d flag to the arg string
        else:
            self.d_ = "-d " + str(d) + " "

            # LATTICE VECTORS OF THE SYSTEM ------------------------

        if (x, y, z) == (None, None, None):  # xyz lattice MUST be supplied.
            raise Exception("x, y and z lattice vectors MUST be given.")
        else:
            self.xyz = "-x " + str(x) + " -y " + str(y) + " -z " + str(z) + " "

        # BOX DIMENSIONS (3, 6, OR 9 FLOATS)
        if box == None:
            self.box_ = ''
        else:
            self.box_ = "-box " + str(box) + " "

        self.pbc_args = self.pbc_ + self.d_ + self.xyz + self.box_
        self.command += self.pbc_args

        self.PBS = True
        print('Periodic setup complete. tags added: ' + self.pbc_args)

    def MembraneSetup(self, l=None, u=None, a=None, au=None, asym=None, hole=None, disc=None, rand=None, bd=None):
        '''collect kwargs to format into command'''

        # LIPID TYPE AND RELATIVE ABUNDANCE ----------------------------
        guide = ["DOPE", "POPG", "DOPG", "CDL1"]

        if l == None:  # you must supply at least one lipid to make the bilayer
            raise Exception(
                "Lipid type must be supplied and match one of the given lipids in lipids.dat or a .itp file inside the current directory. Please fix your -l options.\noptions must be supplied as a list formatted as below:\n\n\t\t\t\toptions=\t['lipidname']\nor for multiple options=\t['lipidname:concentration','lipidname:concenctation',....]\n\nIf you have created a new lipid, you must supply the .itp file in the working directory yourself before running this module, and the file must be named [LIPIDCODE]_martini.itp. for example:\n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"")
        elif type(l) != list:  # it should be a list so we can parse it.
            raise Exception(
                "Lipid type must be supplied and match one of the given lipids in lipids.dat or a .itp file inside the current directory. Please fix your -l options.\noptions must be supplied as a list formatted as below:\n\n\t\t\t\toptions=\t['lipidname']\nor for multiple options=\t['lipidname:concentration','lipidname:concenctation',....]\n\nIf you have created a new lipid, you must supply the .itp file in the working directory yourself before running this module, and the file must be named [LIPIDCODE]_martini.itp. for example:\n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"")

        else:  # if the lipid was supplied and is a list then proceed

            for lipid in l:  # each lipid will be parsed
                if ":" in lipid:  # if theres multiple lipid types, then they will have relative concentrations. we need to separate the lipid code from the numbers.
                    ind = lipid.index(":")
                    lipidcode = lipid[:ind]
                else:
                    lipidcode = lipid

                if lipidcode in guide:  # if we have supplied an itp in the guide, then no additional input is needed.
                    self.l_ = "-l " + lipid + " "

                if lipidcode not in guide:  # if the lipid is a custom lipid then there must be a custom .itp provided.
                    if os.path.isfile(lipidcode + "_martini.itp"):
                        pass

                    elif not os.path.isfile(lipidcode + "_martini.itp"):
                        raise Exception(
                            "Custom Lipids must have associated .itp file supplied inside the working directory and must follow the naming convention:\n\n[LIPIDCODE]_martini.itp. for example:\n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"")
                    else:
                        # if all is well, then do this
                        self.l_ = "-l " + lipid + " "

        if u == None:
            self.u_ = ('')

        elif type(u) != list and type(u) != None:
            raise Exception(
                "Lipid type must be supplied and match one of the given lipids in lipids.dat or a .itp file inside the current directory. Please fix your -u options.\noptions must be supplied as a list formatted as below:\n\n\t\t\t\toptions=\t['lipidname']\nor for multiple options=\t['lipidname:concentration','lipidname:concenctation',....]\n\nIf you have created a new lipid, you must supply the .itp file in the working directory yourself before running this module, and the file must be named [LIPIDCODE]_martini.itp. for example:\n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"")

        else:
            for lipid in u:
                if lipid in guide:
                    self.u_ = "-u " + lipid + " "

                if lipid not in guide:

                    if ":" in lipid:
                        ind = lipid.index(":")
                        lipidcode = lipid[:ind]
                    else:
                        lipidcode = lipid

                    if os.path.isfile(lipidcode + "_martini.itp"):
                        pass

                    elif not os.path.isfile(lipidcode + "_martini.itp"):
                        raise Exception(
                            "Custom Lipids must have associated .itp file supplied inside the working directory and must follow the naming convention:\n\n[LIPIDCODE]_martini.itp. for example:\n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"")
                    else:
                        # if all is well, then do this
                        self.u_ = "-u " + lipid + " "

        # AREA PER LIPID -----------------------------------------------
        if a == None:
            raise Exception("Area per lipid (nm*nm) must be supplied.")
        else:
            self.a_ = "-a " + str(a) + ' '
        if au == None:
            self.au_ = ('')
        else:
            self.au_ = "-au " + str(au) + ' '

        # MEMBRANE SYMMETRY --------------------------------------------
        if asym == None:
            self.asym_ = ''
        else:
            self.asym_ = "-asym " + str(asym) + ' '

        # ARE THERE ANY HOLES????????? ---------------------------------
        if hole == None:
            self.hole_ = ''
        else:
            self.hole_ = "-hole " + str(hole) + ' '

        # MAKE A MEMBRANE DISC -----------------------------------------
        if disc == None:
            self.disc_ = ''
        else:
            self.disc_ = "-disc " + str(disc) + ' '

        # RANDOM KICK SIZE ---------------------------------------------
        if rand == None:
            self.rand_ = ''
        else:
            self.rand_ = "-rand " + str(rand) + ' '

        # BEAD DISTANCE
        if bd == None:
            self.bd_ = ''
        else:
            self.bd_ = "-bd " + str(bd) + ' '

        self.mem_args = self.l_ + self.u_ + self.a_ + self.au_ + self.asym_ + self.hole_ + self.disc_ + self.rand_ + self.bd_

        self.command += self.mem_args
        self.MS = True
        print("Membrane setup complete. arguments added: " + self.mem_args)

    def ProteinSetup(self, center=None, orient=None, rotate=None, od=None, op=None, fudge=None, ring=None, dm=None):
        '''collect kwargs to format into command'''

        # CENTER THE PROTEIN -----------------------------------------
        if center == None:
            self.center_ = ''
        else:
            self.center_ = "-center " + str(center) + ' '

        # ORIENT THE PROTEIN IN THE MEMBRANE -------------------------
        if orient == None:
            self.orient_ = ''
        else:
            self.orient_ = "-orient " + str(orient) + ' '

            # ROTATE THE PROTEIN -----------------------------------------
        if rotate == None:
            self.rotate_ = ''
        else:
            self.rotate_ = "-rotate " + str(rotate) + ' '

        # GRID SPACING FOR DETERMINING THE ORIENTATION ---------------
        if od == None:
            self.od_ = ''
        else:
            self.od_ = "-od " + str(od) + ' '

        # HYDROPHOBIC RATION POWER FOR DETERMINING ORIENTATION -------
        if op == None:
            self.op_ = ''
        else:
            self.op_ = "-op " + str(op) + ' '

        # FUDGE FACTOR
        if fudge == None:
            self.fudge_ = ''
        else:
            self.fudge_ = "-fudge " + str(fudge) + ' '

        # PUT LIPIDS INSIDE THE PROTEIN
        if ring == None:
            self.ring_ = ''
        else:
            self.ring_ = "-ring " + str(ring) + ' '

        # SHIFT THE PROTEIN
        if dm == None:
            self.dm_ = ''
        else:
            self.dm_ = "-dm " + str(dm) + ' '

        self.prot_args = self.center_ + self.orient_ + self.rotate_ + self.od_ + self.op_ + self.fudge_ + self.ring_ + self.dm_

        self.command += self.prot_args

        self.PS = True
        print("Protein Setup complete. arguments added: " + self.prot_args)

    def SolventSetup(self, sol=None, sold=None, solr=None, excl=None, salt=None, charge=None):
        '''collect kwargs to format into command'''

        # SOLVENT TYPE AND RELATIVE ABUNDANCE ---
        if sol == None:
            self.sol_ = ''
        else:
            self.sol_ = "-sol " + str(sol) + ' '

        # SOLVENT DIAMETER -----------------
        if sold == None:
            self.sold_ = ''
        else:
            self.sold_ = "-sold " + str(sold) + ' '

        # SOLVENT RANDOM KICK --------------
        if solr == None:
            self.solr_ = ''
        else:
            self.solr_ = "-solr " + str(solr) + ' '

        # EXCLUSION RANGE --------------
        if excl == None:
            self.excl_ = ''
        else:
            self.excl_ = "-excl " + str(excl) + ' '

        # SALT CONCENTRATION --------------------
        if salt == None:
            self.salt_ = ''
        else:
            self.salt_ = "-salt " + str(salt) + ' '

        # CHARGE OF THE SYSTEM
        if charge == None:
            self.charge_ = ''
        else:
            self.charge_ = "-charge " + str(charge) + ' '

        self.solv_args = self.sol_ + self.sold_ + self.solr_ + self.excl_ + self.salt_ + self.charge_

        self.command += self.solv_args

        self.SlS = True
        print("Solvent Settings complete. arguments added: " + self.solv_args)

    def AdditionalOptions(self, alname=None, alhead=None, allink=None, altail=None):
        'collect kwargs to format into command'

        # ADDITIONAL LIPID NAME
        if alname == None:
            self.alname_ = ''
        else:
            self.alname_ = "-alname " + str(alname) + ' '

        # ADDITIONAL LIPID HEAD SPECIFICATION STRING
        if alhead == None:
            self.alhead_ = ''
        else:
            self.alhead_ = "-alhead " + str(alhead) + ' '

        # ADDITIONAL LIPID LINKER SPECIFICATION STRING
        if allink == None:
            self.allink_ = ''
        else:
            self.allink_ = "-allink " + str(allink) + ' '

        # ADDITIONAL LIPID TAIL SPECIFICATION STRING
        if altail == None:
            self.altail_ = ''
        else:
            self.altail_ = "-altail " + str(altail) + ' '

        self.add_args = self.alname_ + self.alhead_ + self.allink_ + self.altail_

        self.command += self.add_args

        self.AO = True
        print("Additional Options added. arguments added: " + self.add_args)

    def run(self, o):
        '''Run the insane.py command and also change the topol.top file'''

        if self.PBS == False:
            raise Exception("You have not run the periodic setup. Please do so before proceeding.")

        if self.MS == False:
            raise Exception("You have not run the membrane setup. Please do so before proceeding.")

        subprocess.run("python insane.py -o " + str(o) + " -f " + self.system + " -p topol.top " + self.command)
        # os.system("") #might be updated
        print("Insane.py script run successfully. :) ")
        print("You've exceduted the following command:\npython insane.py -o " + str(
            o) + " -f " + self.system + " -p topol.top " + self.command)
        # write the topol.top i have no idea how to do this automatically :( :( :(

        return None

    def Therapy(self):
        '''Just prints a guide to the insanity script'''
        f = open('InsaneHelpTxt.txt', 'r')
        file_contents = f.read()
        print(file_contents)
        f.close()
        return None
