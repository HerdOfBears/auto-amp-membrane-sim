# author jyler menard
# group: Mansbach Lab @ Concordia University, Montreal, QC, Canada

import os
import yaml
import logging

# import .helpers.martinize_custom as martinize
from .helpers import martinize_custom as martinize

__all__ = ["CoarseGrainer"]

class CoarseGrainer():

    def __init__(self):
        # python martinize.py -f GL13Kcoil.pdb  -o gik_1.top -x gik_1.pdb -ff martini22p -ss SSEEEEEEEEESS
        self.configs = {
             "-f": "GL13Kcoil.pdb",
            "-ff": "martini22p",
            # "-ss": "SSEEEEEEEEESS",
            "-ss": ["SSEEEEEEEEESS", "SSEEEEEECCESS"],
            "-out_dir": "output/", # where to write output files
            "-in_dir": "data/"  # where to look for pbd|gro file to coarse-grain
        }
        pass

    def coarse_grain(self, configs:dict):
        """
        Performs coarse-graining of a protein using martinize.py

        Inputs:
            configs: 
                dictionary of options to pass to martinize.py
                Should include include information such as 
                (-f) input protein file (.pdb|.gro), 
                (-o) output topology path (.top), 
                (-x) output cg structure path (.pdb)
                (-ff) the force-field to use
                (-ss) the secondary structure of each residue in the protein.
        """
        logging.info("============== Step 2: Martinize (i.e. coarse-grain) the peptide ================")
        if type(configs) != dict:
            raise TypeError("configs expected to be a dictionary")
        configs = self.output_file_naming(configs=configs) # add output file names to configs
        configs = self.check_configs(configs=configs) # check configs for required options

        # if multiple secondary structures are provided, run martinize for each one
        secondary_structures = configs["-ss"]
        if isinstance(secondary_structures, list):
            for i, second_struct in enumerate(secondary_structures):
                configs[  "-ss"] = second_struct
                configs["-name"] = configs["-f"].split(".")[0] + f"_struct{i}" # add index to input file name
                options = self.run_martinize(configs=configs)
        else:
            configs["-name"] = configs["-f"].split(".")[0] + f"_struct{0}" # add index to input file name
            options = self.run_martinize(configs=configs)

        logging.info("============== Line Break ================")
        return options
    
    def run_martinize(self, configs:dict)->list:

        options = martinize.options # get martinize default options (each option is an object)
        lists   = martinize.lists   # get martinize default, empty lists

        # martinize requires a list of options, not a dictionary. 
        configs = self.convert_config_dict_to_list(configs=configs) # convert configs to list

        logging.info("Martinize parse config file...")
        options = martinize.option_parser(configs, options, lists) # parse config file

        logging.info("Running martinize...")
        martinize.main(options) # run martinize
        logging.info("Martinize complete")

    def load_config_file(self, config_file:str)->dict:
        """ 
        loads a configuration file and returns a dictionary of the options
        """
        with open(config_file, 'r') as f:
            configs = yaml.safe_load(f)
        return configs
    
    def check_configs(self, configs:dict)->dict:
        """
        confirms that the configs dictionary contains required options
        """

        required_options = ["-f", "-o", "-x", "-ff", "-ss"]
        for option in required_options:
            if option not in configs.keys():
                logging.info(f"configs missing required option: {option}")
                raise ValueError(f"configs missing required option: {option}")
            
        logging.info("configs checked successfully")
        return configs
    
    def output_file_naming(self, configs:dict)->dict:
        """
        returns config dictionary with output file names

        Notes:
            relevant martini options:
                (-f) input protein file (.pdb|.gro), 
                (-o) output topology path (.top), 
                (-x) output cg structure path (.pdb)
        """
        input_file = configs["-f"]
        input_file = input_file.split(".")[0] # remove file extension

        output_file_name = input_file + "_martinized"
        configs["-o"] = output_file_name + ".top"
        configs["-x"] = output_file_name + ".pdb"

        logging.warning(f"renamed martini output files to: {output_file_name}")
        return configs

    def convert_config_dict_to_list(self, configs:dict)->list:
        if type(configs) != dict:
            raise TypeError("configs must be a dictionary")
        
        options = []
        for key, value in configs.items():
            options.append(key)
            options.append(value)
        return options
    

if __name__ == "__main__":
    cg = CoarseGrainer()
    
    options = cg.coarse_grain(cg.configs)
    # print(options)