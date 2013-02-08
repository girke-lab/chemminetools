BEGIN;
CREATE TABLE "cms_placeholder" (
    "id" serial NOT NULL PRIMARY KEY,
    "slot" varchar(50) NOT NULL,
    "default_width" smallint CHECK ("default_width" >= 0)
)
;
CREATE TABLE "cms_cmsplugin" (
    "id" serial NOT NULL PRIMARY KEY,
    "placeholder_id" integer REFERENCES "cms_placeholder" ("id") DEFERRABLE INITIALLY DEFERRED,
    "parent_id" integer,
    "position" smallint CHECK ("position" >= 0),
    "language" varchar(15) NOT NULL,
    "plugin_type" varchar(50) NOT NULL,
    "creation_date" timestamp with time zone NOT NULL,
    "changed_date" timestamp with time zone NOT NULL,
    "level" integer CHECK ("level" >= 0) NOT NULL,
    "lft" integer CHECK ("lft" >= 0) NOT NULL,
    "rght" integer CHECK ("rght" >= 0) NOT NULL,
    "tree_id" integer CHECK ("tree_id" >= 0) NOT NULL
)
;
ALTER TABLE "cms_cmsplugin" ADD CONSTRAINT "parent_id_refs_id_e0b32a03" FOREIGN KEY ("parent_id") REFERENCES "cms_cmsplugin" ("id") DEFERRABLE INITIALLY DEFERRED;
CREATE TABLE "cms_page_placeholders" (
    "id" serial NOT NULL PRIMARY KEY,
    "page_id" integer NOT NULL,
    "placeholder_id" integer NOT NULL REFERENCES "cms_placeholder" ("id") DEFERRABLE INITIALLY DEFERRED,
    UNIQUE ("page_id", "placeholder_id")
)
;
CREATE TABLE "cms_page" (
    "id" serial NOT NULL PRIMARY KEY,
    "created_by" varchar(70) NOT NULL,
    "changed_by" varchar(70) NOT NULL,
    "parent_id" integer,
    "creation_date" timestamp with time zone NOT NULL,
    "changed_date" timestamp with time zone NOT NULL,
    "publication_date" timestamp with time zone,
    "publication_end_date" timestamp with time zone,
    "in_navigation" boolean NOT NULL,
    "soft_root" boolean NOT NULL,
    "reverse_id" varchar(40),
    "navigation_extenders" varchar(80),
    "published" boolean NOT NULL,
    "template" varchar(100) NOT NULL,
    "site_id" integer NOT NULL REFERENCES "django_site" ("id") DEFERRABLE INITIALLY DEFERRED,
    "moderator_state" smallint NOT NULL,
    "level" integer CHECK ("level" >= 0) NOT NULL,
    "lft" integer CHECK ("lft" >= 0) NOT NULL,
    "rght" integer CHECK ("rght" >= 0) NOT NULL,
    "tree_id" integer CHECK ("tree_id" >= 0) NOT NULL,
    "login_required" boolean NOT NULL,
    "limit_visibility_in_menu" smallint,
    "publisher_is_draft" boolean NOT NULL,
    "publisher_public_id" integer UNIQUE,
    "publisher_state" smallint NOT NULL
)
;
ALTER TABLE "cms_page_placeholders" ADD CONSTRAINT "page_id_refs_id_b22baae5" FOREIGN KEY ("page_id") REFERENCES "cms_page" ("id") DEFERRABLE INITIALLY DEFERRED;
ALTER TABLE "cms_page" ADD CONSTRAINT "parent_id_refs_id_653a773" FOREIGN KEY ("parent_id") REFERENCES "cms_page" ("id") DEFERRABLE INITIALLY DEFERRED;
ALTER TABLE "cms_page" ADD CONSTRAINT "publisher_public_id_refs_id_653a773" FOREIGN KEY ("publisher_public_id") REFERENCES "cms_page" ("id") DEFERRABLE INITIALLY DEFERRED;
CREATE TABLE "cms_pagemoderator" (
    "id" serial NOT NULL PRIMARY KEY,
    "page_id" integer NOT NULL REFERENCES "cms_page" ("id") DEFERRABLE INITIALLY DEFERRED,
    "user_id" integer NOT NULL REFERENCES "auth_user" ("id") DEFERRABLE INITIALLY DEFERRED,
    "moderate_page" boolean NOT NULL,
    "moderate_children" boolean NOT NULL,
    "moderate_descendants" boolean NOT NULL
)
;
CREATE TABLE "cms_pagemoderatorstate" (
    "id" serial NOT NULL PRIMARY KEY,
    "page_id" integer NOT NULL REFERENCES "cms_page" ("id") DEFERRABLE INITIALLY DEFERRED,
    "user_id" integer REFERENCES "auth_user" ("id") DEFERRABLE INITIALLY DEFERRED,
    "created" timestamp with time zone NOT NULL,
    "action" varchar(3),
    "message" text NOT NULL
)
;
CREATE TABLE "cms_globalpagepermission_sites" (
    "id" serial NOT NULL PRIMARY KEY,
    "globalpagepermission_id" integer NOT NULL,
    "site_id" integer NOT NULL REFERENCES "django_site" ("id") DEFERRABLE INITIALLY DEFERRED,
    UNIQUE ("globalpagepermission_id", "site_id")
)
;
CREATE TABLE "cms_globalpagepermission" (
    "id" serial NOT NULL PRIMARY KEY,
    "user_id" integer REFERENCES "auth_user" ("id") DEFERRABLE INITIALLY DEFERRED,
    "group_id" integer REFERENCES "auth_group" ("id") DEFERRABLE INITIALLY DEFERRED,
    "can_change" boolean NOT NULL,
    "can_add" boolean NOT NULL,
    "can_delete" boolean NOT NULL,
    "can_change_advanced_settings" boolean NOT NULL,
    "can_publish" boolean NOT NULL,
    "can_change_permissions" boolean NOT NULL,
    "can_move_page" boolean NOT NULL,
    "can_moderate" boolean NOT NULL,
    "can_view" boolean NOT NULL,
    "can_recover_page" boolean NOT NULL
)
;
ALTER TABLE "cms_globalpagepermission_sites" ADD CONSTRAINT "globalpagepermission_id_refs_id_2c730e06" FOREIGN KEY ("globalpagepermission_id") REFERENCES "cms_globalpagepermission" ("id") DEFERRABLE INITIALLY DEFERRED;
CREATE TABLE "cms_pagepermission" (
    "id" serial NOT NULL PRIMARY KEY,
    "user_id" integer REFERENCES "auth_user" ("id") DEFERRABLE INITIALLY DEFERRED,
    "group_id" integer REFERENCES "auth_group" ("id") DEFERRABLE INITIALLY DEFERRED,
    "can_change" boolean NOT NULL,
    "can_add" boolean NOT NULL,
    "can_delete" boolean NOT NULL,
    "can_change_advanced_settings" boolean NOT NULL,
    "can_publish" boolean NOT NULL,
    "can_change_permissions" boolean NOT NULL,
    "can_move_page" boolean NOT NULL,
    "can_moderate" boolean NOT NULL,
    "can_view" boolean NOT NULL,
    "grant_on" integer NOT NULL,
    "page_id" integer REFERENCES "cms_page" ("id") DEFERRABLE INITIALLY DEFERRED
)
;
CREATE TABLE "cms_pageuser" (
    "user_ptr_id" integer NOT NULL PRIMARY KEY REFERENCES "auth_user" ("id") DEFERRABLE INITIALLY DEFERRED,
    "created_by_id" integer NOT NULL REFERENCES "auth_user" ("id") DEFERRABLE INITIALLY DEFERRED
)
;
CREATE TABLE "cms_pageusergroup" (
    "group_ptr_id" integer NOT NULL PRIMARY KEY REFERENCES "auth_group" ("id") DEFERRABLE INITIALLY DEFERRED,
    "created_by_id" integer NOT NULL REFERENCES "auth_user" ("id") DEFERRABLE INITIALLY DEFERRED
)
;
CREATE TABLE "cms_title" (
    "id" serial NOT NULL PRIMARY KEY,
    "language" varchar(15) NOT NULL,
    "title" varchar(255) NOT NULL,
    "menu_title" varchar(255),
    "slug" varchar(255) NOT NULL,
    "path" varchar(255) NOT NULL,
    "has_url_overwrite" boolean NOT NULL,
    "application_urls" varchar(200),
    "redirect" varchar(255),
    "meta_description" text,
    "meta_keywords" varchar(255),
    "page_title" varchar(255),
    "page_id" integer NOT NULL REFERENCES "cms_page" ("id") DEFERRABLE INITIALLY DEFERRED,
    "creation_date" timestamp with time zone NOT NULL,
    UNIQUE ("language", "page_id")
)
;
CREATE INDEX "cms_placeholder_slot" ON "cms_placeholder" ("slot");
CREATE INDEX "cms_placeholder_slot_like" ON "cms_placeholder" ("slot" varchar_pattern_ops);
CREATE INDEX "cms_cmsplugin_placeholder_id" ON "cms_cmsplugin" ("placeholder_id");
CREATE INDEX "cms_cmsplugin_parent_id" ON "cms_cmsplugin" ("parent_id");
CREATE INDEX "cms_cmsplugin_language" ON "cms_cmsplugin" ("language");
CREATE INDEX "cms_cmsplugin_language_like" ON "cms_cmsplugin" ("language" varchar_pattern_ops);
CREATE INDEX "cms_cmsplugin_plugin_type" ON "cms_cmsplugin" ("plugin_type");
CREATE INDEX "cms_cmsplugin_plugin_type_like" ON "cms_cmsplugin" ("plugin_type" varchar_pattern_ops);
CREATE INDEX "cms_cmsplugin_level" ON "cms_cmsplugin" ("level");
CREATE INDEX "cms_cmsplugin_lft" ON "cms_cmsplugin" ("lft");
CREATE INDEX "cms_cmsplugin_rght" ON "cms_cmsplugin" ("rght");
CREATE INDEX "cms_cmsplugin_tree_id" ON "cms_cmsplugin" ("tree_id");
CREATE INDEX "cms_page_parent_id" ON "cms_page" ("parent_id");
CREATE INDEX "cms_page_publication_date" ON "cms_page" ("publication_date");
CREATE INDEX "cms_page_publication_end_date" ON "cms_page" ("publication_end_date");
CREATE INDEX "cms_page_in_navigation" ON "cms_page" ("in_navigation");
CREATE INDEX "cms_page_soft_root" ON "cms_page" ("soft_root");
CREATE INDEX "cms_page_reverse_id" ON "cms_page" ("reverse_id");
CREATE INDEX "cms_page_reverse_id_like" ON "cms_page" ("reverse_id" varchar_pattern_ops);
CREATE INDEX "cms_page_navigation_extenders" ON "cms_page" ("navigation_extenders");
CREATE INDEX "cms_page_navigation_extenders_like" ON "cms_page" ("navigation_extenders" varchar_pattern_ops);
CREATE INDEX "cms_page_site_id" ON "cms_page" ("site_id");
CREATE INDEX "cms_page_level" ON "cms_page" ("level");
CREATE INDEX "cms_page_lft" ON "cms_page" ("lft");
CREATE INDEX "cms_page_rght" ON "cms_page" ("rght");
CREATE INDEX "cms_page_tree_id" ON "cms_page" ("tree_id");
CREATE INDEX "cms_page_limit_visibility_in_menu" ON "cms_page" ("limit_visibility_in_menu");
CREATE INDEX "cms_page_publisher_is_draft" ON "cms_page" ("publisher_is_draft");
CREATE INDEX "cms_page_publisher_state" ON "cms_page" ("publisher_state");
CREATE INDEX "cms_pagemoderator_page_id" ON "cms_pagemoderator" ("page_id");
CREATE INDEX "cms_pagemoderator_user_id" ON "cms_pagemoderator" ("user_id");
CREATE INDEX "cms_pagemoderatorstate_page_id" ON "cms_pagemoderatorstate" ("page_id");
CREATE INDEX "cms_pagemoderatorstate_user_id" ON "cms_pagemoderatorstate" ("user_id");
CREATE INDEX "cms_globalpagepermission_user_id" ON "cms_globalpagepermission" ("user_id");
CREATE INDEX "cms_globalpagepermission_group_id" ON "cms_globalpagepermission" ("group_id");
CREATE INDEX "cms_pagepermission_user_id" ON "cms_pagepermission" ("user_id");
CREATE INDEX "cms_pagepermission_group_id" ON "cms_pagepermission" ("group_id");
CREATE INDEX "cms_pagepermission_page_id" ON "cms_pagepermission" ("page_id");
CREATE INDEX "cms_pageuser_created_by_id" ON "cms_pageuser" ("created_by_id");
CREATE INDEX "cms_pageusergroup_created_by_id" ON "cms_pageusergroup" ("created_by_id");
CREATE INDEX "cms_title_language" ON "cms_title" ("language");
CREATE INDEX "cms_title_language_like" ON "cms_title" ("language" varchar_pattern_ops);
CREATE INDEX "cms_title_slug" ON "cms_title" ("slug");
CREATE INDEX "cms_title_slug_like" ON "cms_title" ("slug" varchar_pattern_ops);
CREATE INDEX "cms_title_path" ON "cms_title" ("path");
CREATE INDEX "cms_title_path_like" ON "cms_title" ("path" varchar_pattern_ops);
CREATE INDEX "cms_title_has_url_overwrite" ON "cms_title" ("has_url_overwrite");
CREATE INDEX "cms_title_application_urls" ON "cms_title" ("application_urls");
CREATE INDEX "cms_title_application_urls_like" ON "cms_title" ("application_urls" varchar_pattern_ops);
CREATE INDEX "cms_title_page_id" ON "cms_title" ("page_id");
COMMIT;
