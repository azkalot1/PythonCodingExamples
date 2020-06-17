class ScikitData:
    # just a dummy example
    def __init__(self):
        super().__init__()

    def pair_plot(self, vars=range(3), hue=None):
        return pairplot(pd.DataFrame(self.data), vars=vars, hue=hue, kind="reg")
    
    @classmethod
    # allows to do 
    # dataset_generator = ScikitData.get_generator(["diabetes", "iris"])
    def get_generator(cls, dataset_names):
        return map(cls, dataset_names)



