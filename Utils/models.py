from django.core.exceptions import ValidationError
from django.db import models


KNOWN_SPECIAL_FIELDS_TO_IGNORE = ['id', 'time_created_on_server', 'time_last_updated_on_server', 'preferred_order']


class SortableModel(models.Model):
    preferred_order = models.PositiveSmallIntegerField(default=32767)

    class Meta:
        abstract = True
        ordering = ['preferred_order']


class BigIdModel(models.Model):
    id = models.BigAutoField(primary_key=True)

    class Meta:
        abstract = True


class LogCreationTimeModel(models.Model):
    time_created_on_server = models.DateTimeField(auto_now_add=True)  # TODO: Make this not editable in the admin (editable=False makes it disappear)

    class Meta:
        abstract = True


class LogLastUpdatedTimeModel(models.Model):
    time_last_updated_on_server = models.DateTimeField(auto_now=True)  # TODO: Make this not editable in the admin (editable=False makes it disappear)

    class Meta:
        abstract = True


class CustomDeleteQuerySet(models.QuerySet):
    def delete(self):
        for obj in self:
            obj.delete()


class CustomDeleteManager(models.Manager):
    """
    Any model that overrides the delete method should use this Manager
    objects = CustomDeleteManager()
    """
    def get_queryset(self):
        queryset = CustomDeleteQuerySet(self.model, using=self._db)
        return queryset


class CustomSaveQuerySet(models.QuerySet):
    def create(self, **kwargs):
        raise NotImplementedError("You must call the model save on this model!")

    def bulk_create(self, objs, batch_size=None):
        raise NotImplementedError("You must call the model save on this model!")

    def get_or_create(self, defaults=None, **kwargs):
        raise NotImplementedError("You must call the model save on this model!")

    def update_or_create(self, defaults=None, **kwargs):
        raise NotImplementedError("You must call the model save on this model!")

    def update(self, force=False, **kwargs):
        if force:
            return super(CustomSaveQuerySet, self).update(**kwargs)
        else:
            raise NotImplementedError("You must call the model save on this model!")


class CustomSaveManager(models.Manager):
    """
    Any model that overrides the save method should use this Manager
    objects = CustomSaveManager()
    """
    def get_queryset(self):
        queryset = CustomSaveQuerySet(self.model, using=self._db)
        return queryset


class CustomSaveAndDeleteQuerySet(CustomSaveQuerySet, CustomDeleteQuerySet):
    pass


class CustomSaveAndDeleteManager(models.Manager):
    """
    Any model that overrides the save and delete methods should use this Manager
    objects = CustomSaveManager()
    """
    def get_queryset(self):
        queryset = CustomSaveAndDeleteQuerySet(self.model, using=self._db)
        return queryset


class UnEditableModel(models.Model):
    objects = CustomSaveManager()

    def save(self, *args, **kwargs):
        if self.id:
            raise ValidationError("You cannot edit instances of %s!" % self._meta.object_name)
        return super(UnEditableModel, self).save(*args, **kwargs)

    class Meta:
        abstract = True
