cd /srv/www
rsync -pgrtuzzvl -e ssh --progress htdocs musil:/srv/www
cd -


