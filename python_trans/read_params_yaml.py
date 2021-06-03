import yaml
class read_params_yaml:
    def __init__(self, file_name):
        
        with open(file_name) as file:
            yaml_dic = yaml.load(file, Loader=yaml.FullLoader)
        def sample_token(self,**response):
            for k,v in response.items():
                if isinstance(v,dict):
                    self.__dict__[k] = sample_token(self,**v)
                else:
                     self.__dict__[k] = v
        sample_token(self,**yaml_dic)