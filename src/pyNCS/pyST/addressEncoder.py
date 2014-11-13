import numpy as np

class AddressEncoder:
    def __init__(self, addr_conf, addr_spec, nBits, nBitsTotal, addr_pinconf):
        '''
        This class carries the functions to translate between physical and human readable addresses.
        '''
        #Construct Extract, Construct
        #self.fc, fc_field = self._stas_create_construct(addr_conf, addr_pinconf)
        self.fe, fe_field = self._stas_create_extract(addr_conf, addr_pinconf)

        nDims = len(addr_conf)

        self.fc_field_dict = {}
        self.fe_field_dict = {}

        id_list = extract_id_list(addr_pinconf) # List of human readable fields
        #for i in range(len(addr_pinconf)):
        #    self.fc_field_dict[id_list[i]] = fc_field[i]

        id_list = extract_id_list(addr_conf)
        for i in range(len(addr_conf)):
            self.fe_field_dict[id_list[i]] = fe_field[i]

    def __getitem__(self, item):
        if item in self.fc_field_dict.keys():
            return self.fc_field_dict[item]
        elif item in self.fe_field_dict.keys():
            return self.fe_field_dict[item]
        else:
            raise AttributeError("There is no such field or pin %s" % item)

    def encode(self, addr):
        return np.vstack(self.fc(addr))

    def decode(self, addr):
        return np.vstack(self.fe(addr))

    def encode_field(self, field, addr):
        return np.array(self[field](*addr), 'uint32')

    def decode_field(self, field, addr):
        return np.array(self[field](*addr), 'uint32')


    def _stas_create_extract(self, addr_conf, addr_pinconf):
        """
        Parsing the Extraction functions
        of the type Physical -> Human
        """
        #TODO: These funcitons are purely dependent on the order in addrConf
        arg_str = ''.join([i + ', ' for i in extract_id_list(addr_pinconf)])
        fe_field = [eval('lambda ' + arg_str + ':' + s['f']) for s in addr_conf]
        fe = eval('lambda ' + arg_str + ':' + "[" + ''.join([s['f'] +
            ', ' for s in addr_conf]) + ']')
        return lambda x: fe(*x), fe_field
    
    
    def _stas_create_construct(self, addr_conf, addr_pinconf):
        """
        Parsing the Construction functions
        of the type Human -> Physical
        """
        #TODO: These funcitons are purely dependent on the order in addrConf
        arg_str = ''.join([i + ', ' for i in extract_id_list(addr_conf)])
        fc_field = [eval('lambda ' + arg_str + ':' + v['f']) for v in addr_pinconf]
        fc = eval('lambda ' + arg_str + ':' + "[" + ''.join([s['f'] +
            ', ' for s in addr_pinconf]) + ']')
        return lambda x: fc(*x), fc_field

def extract_id_list(addr_conf):
    id_list = []
    for i in addr_conf:
        id_list.append(i['id'])
    return id_list
