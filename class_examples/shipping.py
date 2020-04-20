class ShippingContainer:
    """
    Simple example of class
    """
    next_serial = 1337  # Class attribute
    @staticmethod
    def _get_next_serial_staticmethod():
        result = ShippingContainer.next_serial
        ShippingContainer.next_serial += 1
        return result
        # Assosiate _get_next_serial with a class, not instance
        # of a class

    @classmethod
    def _get_next_serial_classmethod(cls):
        result = cls.next_serial
        cls.next_serial += 1
        return result
        # Assosiate _get_next_serial with a class, not instance
        # of a class

    @classmethod
    def create_empty(cls, owner_code, *args, **kwargs):
        # Usefull for derived classes
        return cls(owner_code, content=None, *args, **kwargs)
        # Will call class and makes an instance
        # Can be used as ShippingContainer.create_empty("YML")

    @staticmethod
    def _make_bic_code(owner_code, serial):
        return iso6346.create(
            owner_code=owner_code,
            serial=str(serial).zfill(6))
        # Let pretend it was imported

    def __init__(self, owner_code, content):
        self.owner_code = owner_code  # It will be an instance attribute
        self.content = content  # Another instance attribute
        self.serial = ShippingContainer.next_serial  # To adress class atribute
        self.bic = self._make_bic_code(
            owner_code=owner_code,
            serial=ShippingContainer._get_next_serial()
        )
        # we can write ShippingContainer._make_bic_code, however, there
        # will be no override from child class when called
        ShippingContainer.next_serial += 1
        # Alternatively self.serial = self._get_next_serial()
        # Will increment serial number as we make new objects of this class


# Now add a child class

class RefrigeratedShippingContainer(ShippingContainer):

    MAX_CELSIUS = 4.0
    # Will apply to all instances of RefrigeratedShippingContainer
    @staticmethod
    def _c_to_f(celsius):
        return celsius * 9/5 + 32
 
    @staticmethod
    def _f_to_c(fahrenheit):
        return(fahrenheit - 32) * 5/9

    @staticmethod
    def _make_bic_code(owner_code, serial):
        return iso6346.create(
            owner_code=owner_code,
            category='R',
            serial=str(serial).zfill(6))
    # This overwrites parent staticmethod from ShippingContainer
    # Class methods will behave polymorphicaly
    # so RefrigeratedShippingContainer.create_empty() will
    # make an instance of RefrigeratedShippingContainer
    # because it will pass appropriate class

    def __init__(self, owner_code, content, celsius):
        super().__init__(owner_code, content)
        # Will call init of parent class
        self.celsius = celsius

    @property
    def celsius(self):
        return self._celsius
        # Now we can't modify by c.celsius directly!

    # node let define setter
    @celsius.setter
    def celsius(self, value):
        if value > RefrigeratedShippingContainer.MAX_CELSIUS:
            # We can't just call MAX_CELSIUS, it is outside of the scope
            raise ValueError("Too hot")       
        self._celsius = value
    # Example of usage
    # r = RefrigeratedShippingContainer('YML','prawn', celsius=-10)
    # r.celsius = -5 <- Will call validation method again
    # Which will allow proper setting

    # Now add property for Fahrenheit

    @property
    def fahrenheit(self):
        return RefrigeratedShippingContainer._c_to_f(self.celsius)
    # This will set a new property, fahrenheit in init???

    @fahrenheit.setter
    def fahrenheit(self, value):
        self.celsius = RefrigeratedShippingContainer._c_to_f(self.celsius)
        # it invokes @celsius.setter 
    # If we set the temperatire using r.fahrenheit
    # it will invokes fahrenheit setter and update celsius as well
    
    # overriding the properties also works


class HeatedRefrigeratedShippingContainer(RefrigeratedShippingContainer):

    MIN_CELSIUS = -20

    # Now to set celsius we need to override celcius setter
    # However we must explicitly set which setter we override
    # However we don't need to override the property itself
    @RefrigeratedShippingContainer.celsius.setter
    def celsius(self, value):
        if not (HeatedRefrigeratedShippingContainer.MIN_CELSIUS
                <= value
                <= RefrigeratedShippingContainer.MAX_CELSIUS):
            raise ValueError('Out of range')
        self._celsius = value
