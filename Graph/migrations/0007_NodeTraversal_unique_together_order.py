# Generated by Django 2.2.1 on 2019-09-03 16:07

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('Graph', '0006_Path_zoomlevel'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='nodetraversal',
            unique_together={('path', 'order')},
        ),
    ]
