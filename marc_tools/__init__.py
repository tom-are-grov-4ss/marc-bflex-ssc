import pkg_resources
import logging

# from . import general
# from . import scf


__version__ = pkg_resources.get_distribution('marc_tools').version

logging.basicConfig(filename=r'./Marc_tools.log', filemode='w',
                    level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
