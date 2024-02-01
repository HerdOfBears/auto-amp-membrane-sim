# the class needs to take all the possible options for the insane script, then use the os package to run the code
import os
import subprocess
import pandas as pd
import logging

__all__ = ["Insanity"]


class Insanity():
    '''Wrapper function to take options from the user to run the insane script. it is possible to import the insane
    package but they don't have any documentation so using it is impossible.

    Heed the instructions!
    '''

    def __init__(self, sys, Defaults=True):
        self.nochains = 20 #this is hard coded rn because im not sure how to add this automatically from the prev step
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
        self.u_ = ''
        self.l_ = ''
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
        self.include = ''
        # check for insane.py in the folder
        if os.path.isfile('insane.py'):
            pass
        else:
            raise Exception("Missing script file: insane.py. Please locate and insert into directory")

        logging.info('Initialized insane script object. If you need help, Please run the therapy method!')

    def add_separate_chains_to_file(self, filename, nochains, outname):
        ''' written by Re Mansbach

        Adds the chains previously added to the structure to the top file for indexing reasons. this is necessary IFF the user wants to be able to index their chains easier after the simulation step. only run this AFTER the self.run function was completed and your topol.top default file has been created.

        available:

        filename (str):
            the original top file created after the insane step (should be called topol.top, but can be changed manually)

        nochains (float):
            number of chains added to the membrane structure, taken from previous step (soon to be implemented)

        outname (str):
            output filename (can overwrite topol.top if you want, but recommended to save as something unique)

        returns

        outputs new .top file with added indexing for analysis ease, and include statements added after the proper martini include statement. check working directory for file.
        '''
        letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789"
        if nochains > len(letters):
            raise ValueError("too many chains")
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

        outf = open(outname, 'w')
        for line in lines:
            if "#include" in line:
                for incl in self.include:
                    outf.write(incl + "\n")
            elif "Protein" in line and "INSANE" not in line:
                for i in range(nochains):
                    outf.write("Protein {}\t\t 1\n".format(letters[i]))
            else:
                outf.write(line)
        outf.close()

    def periodic_boundary_setup(self, pbc=None, d=None, dz=None, x=None, y=None, z=None, box=None, n=None):
        """collect kwargs and format into command to satisfy the the insane.py script command

        available:

        pbc (string):
            periodic boundary conditions
            must be supplied, and be one of the given options: ['hexagonal', 'rectangular', 'square', 'cubic', 'optimal', 'keep'] otherwise this will raise exception.

        d (string):
            distance between periodic images

        x, y, z (string, float, int):
            lattice vectors of the system
            must be supplied, otherwise this will raise an exception.

        returns:

        (technically this returns None, but it will enable you to access self.pbc_args which contains the flags inputted for this section of the protocol.)
        """
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
        logging.info('Periodic setup complete. tags added: ' + self.pbc_args)

    def membrane_setup(self, l=None, u=None, a=None, au=None, asym=None, hole=None, disc=None, rand=None, bd=None):
        """collect kwargs and format into command to satisfy the insane.py script command

        available:

        l, u (list):
            Lipid type and relative abundance (NAME[:#]). Must be formatted correctly, and corresponding .itp file must be supplied for user-made lipids. (-u refers to the upper leaflet iff the tag was supplied. otherwise, the -l flag will duplicate the upper and lower leaflet)

        a, au (float):
            Area per lipid (nm*nm) (au refers to the upper leaflet iff -u tag was supplied.)

        asym (float):
            Membrane asymmetry (number of lipids)

        hole (float):
            Make a hole in the membrane with specified (given) radius

        disc (float):
            Make a membrane disc with specified radius

        rand (float):
            Random kick size (maximum atom displacement)

        bd (float):
            Bead distance unit for scaling z-coordinates (nm)

        returns:

        (technically this returns None, but it will enable you to access self.mem_args which contains the flags inputted for this section of the protocol.)
        """

        # set up a list to collect the #include lines

        self.include = ['#include \"martini.itp\"']

        # LIPID TYPE AND RELATIVE ABUNDANCE ----------------------------

        if l == None:  # you must supply at least one lipid to make the bilayer
            raise Exception(
                "Lipid type must be supplied and match one of the given lipids in all_lipid_guide.dat or a .itp file inside the current directory. Please fix your -l options.\noptions must be supplied as a list formatted as below:\n\n\t\t\t\toptions=\t['lipidname']\nor for multiple options=\t['lipidname:concentration','lipidname:concenctation',....]\n\nIf you have created a new lipid, you must supply the .itp file in the working directory yourself before running this module, and the file must be named [LIPIDCODE]_martini.itp. for example:\n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"")
        elif type(l) != list:  # it should be a list so we can parse it.
            raise Exception(
                "Lipid type must be supplied and match one of the given lipids in all_lipid_guide.dat or a .itp file inside the current directory. Please fix your -l options.\noptions must be supplied as a list formatted as below:\n\n\t\t\t\toptions=\t['lipidname']\nor for multiple options=\t['lipidname:concentration','lipidname:concenctation',....]\n\nIf you have created a new lipid, you must supply the .itp file in the working directory yourself before running this module, and the file must be named [LIPIDCODE]_martini.itp. for example:\n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"")

        else:  # if the lipid was supplied and is a list then proceed
            lipid_guide = pd.read_csv('all_lipids_v2.0/all_lipid_guide.csv')

        for lipid in l:  # each lipid will be parsed
            if ":" in lipid:  # if theres multiple lipid types, then they will have relative concentrations. we need to separate the lipid code from the numbers.
                ind = lipid.index(":")
                lipidcode = lipid[:ind]
            else:
                lipidcode = lipid
            # check if the itp is in the filefolder first.
            if lipidcode not in list(lipid_guide[
                                         'Lipid_name']):  # if the lipid is a custom lipid then there must be a custom .itp provided.
                if os.path.isfile(lipidcode + "_martini.itp"):
                    logging.warning(
                        "File found:" + lipidcode + "_martini.itp\nThis lipid is not implemented in our package, you will use it at your own discretion")
                    self.l_ += "-l " + lipid + " "
                    self.include.append("#include \"" + lipidcode + "_martini.itp\"")

                elif not os.path.isfile(lipidcode + "_martini.itp"):
                    raise Exception('''Custom Lipids must have associated .itp file supplied inside the working 
                    directory and must follow the naming convention:\n\n[LIPIDCODE]_martini.itp. for example:
                    \n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"\n\nThe offending peptide is: {}'''.format(
                        lipidcode))

            elif lipidcode in list(lipid_guide[
                                       'Lipid_name']):  # if we have supplied an itp in the guide, then no additional input is needed.
                # check for user provided itp file here
                self.l_ += "-l " + lipid + " "
                if os.path.isfile(lipidcode + "_martini.itp"):
                    # use that itp instead of the one in the guide.
                    logging.warning(
                        "File found:" + lipidcode + "_martini.itp\nWarning: we will use this file instead of the one present in our guide.")
                    self.include.append("#include \"" + lipidcode + "_martini.itp\"")
                elif not os.path.isfile(lipidcode + "_martini.itp"):
                    logging.info("Using default itp files.")
                    self.include.append("#include \"" + lipidcode + "_martini.itp\"")

        if u == None:
            self.u_ = ''

        else:  # if the lipid was supplied and is a list then proceed
            lipid_guide = pd.read_csv('all_lipids_v2.0/all_lipid_guide.csv')

            for lipid in u:  # each lipid will be parsed
                if ":" in lipid:  # if theres multiple lipid types, then they will have relative concentrations. we need to separate the lipid code from the numbers.
                    ind = lipid.index(":")
                    lipidcode = lipid[:ind]
                else:
                    lipidcode = lipid
                # check if the itp is in the filefolder first.

                if lipidcode not in list(lipid_guide[
                                             'Lipid_name']):  # if the lipid is a custom lipid then there must be a custom .itp provided.
                    if os.path.isfile(lipidcode + "_martini.itp"):
                        logging.warning(
                            "File found:" + lipidcode + "_martini.itp\nThis lipid is not implemented in our package, you will use it at your own discretion")
                        self.u_ += "-u " + lipid + " "
                        self.include.append("#include \"" + lipidcode + "_martini.itp\"")

                    elif not os.path.isfile(lipidcode + "_martini.itp"):
                        raise Exception('''Custom Lipids must have associated .itp file supplied inside the working 
                        directory and must follow the naming convention:\n\n[LIPIDCODE]_martini.itp. for example:
                        \n\nLIPIDCODE\t\tFILENAME\nABCD\t\t\"ABCD_martini.itp\"''')

                elif lipidcode in list(lipid_guide[
                                           'Lipid_name']):  # if we have supplied an itp in the guide, then no additional input is needed.
                    # check for user provided itp file here
                    self.u_ += "-u " + lipid + " "
                    if os.path.isfile(lipidcode + "_martini.itp"):
                        # use that itp instead of the one in the guide.
                        logging.warning(
                            "File found:" + lipidcode + "_martini.itp\nWarning: we will use this file instead of the one present in our guide.")
                        self.include.append("#include \"" + lipidcode + "_martini.itp\"")
                    elif not os.path.isfile(lipidcode + "_martini.itp"):
                        logging.info("Using default itp files.")
                        self.include.append("#include \"" + lipidcode + "_martini.itp\"")

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
        logging.info("Membrane setup complete. arguments added: " + self.mem_args)

    def protein_setup(self, center=None, orient=None, rotate=None, od=None, op=None, fudge=None, ring=None, dm=None):
        '''collect kwargs and format into command to satisfy the insane.py script command

        available:

        center (bool):
            Center the protein on z

        orient (bool):
            Orient protein in membrane

        rotate (float):
            Rotate protein (random|princ|angle(float))

        od (float):
            Grid spacing for determining orientation

        -op (float):
            Hydrophobic ratio power for determining orientation

        fudge (float):
            Fudge factor for allowing lipid-protein overlap

        ring (bool):
            Put lipids inside the protein

        dm (float):
            Shift protein with respect to membrane

        returns:

        (technically this returns None, but it will enable you to access self.prot_args which contains the flags inputted for this section of the protocol.)

       '''

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
        logging.info("Protein Setup complete. arguments added: " + self.prot_args)

    def solvent_setup(self, sol=None, sold=None, solr=None, excl=None, salt=None, charge=None):
        '''collect kwargs and format into command to satisfy the insane.py script command

        available:

        sol (NAME[:float]):
            Solvent type and relative abundance (NAME[:#])

        sold (float):
            Solvent diameter

        solr (float):
            Solvent random kick

        excl (float):
            Exclusion range (nm) for solvent addition relative to membrane center

        salt(float):
            Salt concentration

        charge(bool):
            Charge of system. Set to auto to infer from residue names

        returns:

        (technically this returns None, but it will enable you to access self.solv_args which contains the flags inputted for this section of the protocol.)

       '''

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
        logging.info("Solvent Settings complete. arguments added: " + self.solv_args)

    def additional_options(self, alname=None, alhead=None, allink=None, altail=None):
        '''collect kwargs and format into command to satisfy the insane.py script command. This command contains additional lipids. see insane script documentation for more information.

        returns:
        (technically this returns None, but it will enable you to access self.add_args which contains the flags inputted for this section of the protocol.)

       '''

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
        logging.info("Additional Options added. arguments added: " + self.add_args)

    def run(self, output_filename):
        '''Run the insane.py command to create the membrane system and also change the topol.top file to include all relevant information.'''

        if self.PBS == False:
            raise Exception("You have not run the periodic setup. Please do so before proceeding.")

        if self.MS == False:
            raise Exception("You have not run the membrane setup. Please do so before proceeding.")

        subprocess.run(
            "python insane.py -o " + str(output_filename) + " -f " + self.system + " -p topol.top " + self.command)
        logging.info("Insane.py script run successfully. :) ")
        logging.info("You've excecuted the following command:\npython insane.py -o " + str(
            output_filename) + " -f " + self.system + " -p topol.top " + self.command)
        self.add_separate_chains_to_file('topol.top',self.nochains,'topol_modified.top')

        return None

    def therapy(self):
        '''Just prints a guide to the insanity script'''
        f = open('InsaneHelpTxt.txt', 'r')
        file_contents = f.read()
        print(file_contents)
        f.close()
        return None
